# datasources.py - data sources handling for the cancermuts package
# (c) 2018 Matteo Tiberti <matteo.tiberti@gmail.com>
# This file is part of cancermuts
#
# cancermuts is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Nome-Programma is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Nome-Programma.  If not, see <http://www.gnu.org/licenses/>.


"""
datasources classes --- :mod:`cancermuts.datasources`
================================================================
Classes to interrogate data sources and annotate various data

"""

import time
import requests as rq
from bioservices.uniprot import UniProt as bsUniProt
from .core import Sequence, Mutation
from .properties import *
from .metadata import *
from .log import *
import re
import sys
import os
import csv
import myvariant
import pyliftover

reload(sys)
sys.setdefaultencoding('utf8')

class Source(object):
    def __init__(self, name, version, description):
        self.name = name
        self.version = version
        self.description = description
    def get_sequence(self, gene_id):
        return None
    def get_mutations(self, ene_id):
        return None
    def get_position_properties(self, position):
        return None
    def get_mutation_properties(self, mutation):
        return None

class DynamicSource(Source, object):
    def __init__(self, *args, **kwargs):
        super(DynamicSource, self).__init__(*args, **kwargs)

class StaticSource(Source, object):
    def __init__(self, *args, **kwargs):
        super(StaticSource, self).__init__(*args, **kwargs)

class UniProt(DynamicSource, object):
    @logger_init
    def __init__(self, *args, **kwargs):
        description = "Uniprot knowledge-base"
        super(UniProt, self).__init__(name='UniProt', version='1.0', description=description)
        self._uniprot_service = bsUniProt()

    def get_sequence(self, gene_id, upid=None):
        if upid is None:
            self.log.info("retrieving UniProt ID for human gene %s" % gene_id)
            upids = self._uniprot_service.search("gene:%s AND organism:human" % gene_id, columns="entry name").split()[2:]
            upid = upids[0]
            if len(upids) > 1:
                self.log.warning("the following UniProt entries were found for gene %s: %s; will use %s" %(gene_id, ', '.join(upids), upid))
            else:
                self.log.info("will use Uniprot ID %s" % upid)
            aliases = {}
        else:
            aliases = {'uniprot':upid}

        self.log.info("retrieving sequence for UniProt sequence for Uniprot ID %s" % gene_id)

        try:
            sequence = self._uniprot_service.get_fasta_sequence(str(upid))
        except:
            self.log.error("failed retrieving sequence for Uniprot ID %s" % upid)
            return None



        return Sequence(gene_id, sequence, self, aliases=aliases)

class cBioPortal(DynamicSource, object):
    @logger_init
    def __init__(self):
        description = "cBioPortal"
        super(cBioPortal, self).__init__(name='cBioPortal', version='1.0', description=description)
        self._requests_url = 'http://www.cbioportal.org/webservice.do'
        self._cache_cancer_studies = None
        self._cache_genetic_profiles = None
        self._cache_case_sets = None
        self._cache_profile_data = None
        self._cancer_types = {}
        self._get_cancer_types()

    def add_mutations(self, sequence, cancer_studies=None, metadata=[], mainisoform=1, mapisoform=None):
        mutations, out_metadata = self._get_profile_data(sequence.gene_id, cancer_studies=cancer_studies, metadata=metadata, mainisoform=mainisoform, mapisoform=mapisoform)
        unique_mutations = list(set(mutations))
        self.log.info("unique mutations found in %s: %s" % (self.name, ", ".join(sorted(unique_mutations, key=lambda x: int(x[1:-1])))))
        for m_idx,m in enumerate(unique_mutations):
            wt = m[0]
            mut = m[-1]
            num = int(m[1:-1])
            if wt == mut:
                self.log.info("synonymous mutation %s discarded" % m)
                continue

            mutation_indices = [i for i,x in enumerate(mutations) if x == m]

            try:
                site_seq_idx = sequence.seq2index(num)
            except:
                self.log.warning("mutation %s is outside the protein sequence; it will be skipped")
                continue

            position = sequence.positions[site_seq_idx]

            if position.wt_residue_type != wt:
                self.log.warning("for mutation %s, residue %d is %s in wild-type sequence; it will be skipped" %(m, wt, position.wt_residue_type))
                continue

            mutation_obj = Mutation(sequence.positions[site_seq_idx],
                                    mut,
                                    [self])
            for md in metadata:
                mutation_obj.metadata[md] = []
                for mi in mutation_indices:
                    if out_metadata[md][mi] is not None:
                        tmp_md = [self] + out_metadata[md][mi]
                        this_md = metadata_classes[md](*tmp_md)
                        mutation_obj.metadata[md].append(this_md)

            self.log.debug("adding mutation %s" % str(mutation_obj))
            sequence.positions[site_seq_idx].add_mutation(mutation_obj)

    def _get_cancer_types(self):
        self.log.info("fetching cancer types")
        try:
            response = rq.post(self._requests_url, {'cmd':'getTypesOfCancer'}).text
        except:
            self.log.error("failed retrieving cancer types")
            self._cancer_types = None
            return    
        tmp = response.split("\n")
        data = []
        for line in tmp[:-1]:
            data.append(tuple(map(str, line.strip().split("\t"))))
        self._cancer_types = dict(data)

    def _get_cancer_studies(self):
        self.log.info("retrieving cancer studies...")
        self._cache_cancer_studies = []
        try:
            response = rq.post(self._requests_url, {'cmd':'getCancerStudies'}).text
        except:
            self.log.error("failed retrieving cancer studies")
            self._cache_cancer_studies = None
            return
        tmp = response.split("\n")
        for line in tmp[1:]:
            if line:
                tmp2 = line.split("\t")
                self._cache_cancer_studies.append(tuple(tmp2))

    def _get_genetic_profiles(self, cancer_studies=None):
        self.log.info("retrieving genetic profiles")
        if cancer_studies is not None:
            self.log.info("user-provided cancer studies will be used")
            self._cache_cancer_studies = cancer_studies
        else:
            self.log.info("all available cancer studies will be used")
            self._get_cancer_studies()

        self._cache_genetic_profiles = {}
        for cancer_study_id in zip(*self._cache_cancer_studies)[0]:
            self._cache_genetic_profiles[cancer_study_id] = []
            self.log.debug("fetching genetic profiles for cancer study %s" % cancer_study_id)
            try:
                response = rq.post(self._requests_url, {'cmd':'getGeneticProfiles',
                                                        'cancer_study_id':cancer_study_id}).text
            except:
                self.log.error("failed to fetch genetic profiles for study %s" % cancer_study_id)
                continue
            tmp = response.split("\n")
            for line in tmp[1:]:
                if line:
                    tmp2 = line.split("\t")
                    if tmp2[4] == 'MUTATION' or tmp2[4] == 'MUTATION_EXTENDED':
                        self._cache_genetic_profiles[cancer_study_id].append((tmp2[0],tmp2[4]))

    def _get_case_sets(self, cancer_studies=None):
        self.log.info("retrieving case sets")
        if cancer_studies is not None:
            self.log.info("user-provided cancer studies will be used")
            self._cache_cancer_studies = cancer_studies
        else:
            self.log.info("all available cancer studies will be used")
            self._get_cancer_studies()

        self._cache_case_sets = {}
        for cancer_study_id in zip(*self._cache_cancer_studies)[0]:
            self._cache_case_sets[cancer_study_id] = []
            self.log.debug("fetching case sets for cancer study %s" % cancer_study_id)
            try:
                response = rq.post(self._requests_url, {'cmd':'getCaseLists',
                                                        'cancer_study_id':cancer_study_id}).text
            except:
                self.log.error("failed to fetch case sets for study %s" % cancer_study_id)
                continue

            tmp = response.split("\n")
            for line in tmp[1:]:
                if line:
                    tmp2 = line.split("\t")
                    self._cache_case_sets[cancer_study_id].append(tmp2[0])

    def _get_profile_data(self, gene_id, cancer_studies=None, metadata=[], mapisoform=None, mainisoform=1):

        self.log.info("fetching profile data")
        mut_regexp = '[A-Z][0-9]+[A-Z]$'
        mut_prog = re.compile(mut_regexp)
        split_regexp = '\s+|,'
        split_prog = re.compile(split_regexp)

        out_metadata = dict(zip(metadata, [list() for i in range(len(metadata))]))

        self._get_genetic_profiles(cancer_studies=cancer_studies)
        self._get_case_sets(cancer_studies=cancer_studies)

        if len(self._cache_cancer_studies[0]) == 3:
            for i,c in enumerate(self._cache_cancer_studies):
                self._cache_cancer_studies[i] = c + tuple([mainisoform])    

        responses = []
        mutations = []
        do_cancer_type = False
        if 'cancer_type' in metadata:
            out_metadata['cancer_type'] = []
            do_cancer_type = True
        do_cancer_study = False
        if 'cancer_study' in metadata:
            out_metadata['cancer_study'] = []
            do_cancer_study = True
        do_genomic_coordinates = False
        if 'genomic_coordinates' in metadata:
            out_metadata['genomic_coordinates'] = []
            do_genomic_coordinates = True
        do_genomic_mutations = False
        if 'genomic_mutations' in metadata:
            out_metadata['genomic_mutations'] = []
            do_genomic_mutations = True

        for cancer_study_id, short_desc, long_desc, cancer_study_isoform in self._cache_cancer_studies:
            self.log.debug("fetching profile data for study %s" % cancer_study_id)
            if do_cancer_type:
                cancer_study_suffix = cancer_study_id.split("_")[0]
                try:
                    cancer_type = self._cancer_types[cancer_study_suffix]
                except KeyError:
                    self.log.error("cancer type %s not found - cancer study will be used instead")
                    cancer_type = cancer_study_suffix
            cancer_study_mutation_ids =     [ x[0] for x in self._cache_genetic_profiles[cancer_study_id] if x[1] == 'MUTATION' ]
            cancer_study_mutation_ext_ids = [ x[0] for x in self._cache_genetic_profiles[cancer_study_id] if x[1] == 'MUTATION_EXTENDED' ]

            for case_set in self._cache_case_sets[cancer_study_id]:
                self.log.debug("doing case set %s for genetic profiles %s" % (case_set, ",".join(cancer_study_mutation_ids)))
                if len(cancer_study_mutation_ids) > 0:
                    try:
                        response = rq.post(self._requests_url, {'cmd':'getProfileData',
                                                                'case_set_id':case_set,
                                                                'genetic_profile_id':",".join(cancer_study_mutation_ids),
                                                                'gene_list':gene_id}).text
                    except:
                        log.error("failed to fetch profile data for these case set and genetic profiles")
                        return None

                    for line in response.strip().split("\n"):
                        if not line.startswith("#") and not line.startswith("GENE_ID"):
                            tmp = split_prog.split(line)[2:]
                            tmp2 = filter(lambda x: x != 'NaN', tmp)
                            tmp3 = filter(lambda x: mut_prog.match(x), tmp2)
                            fixed_isoform = [self._convert_isoform(m, mainisoform, cancer_study_isoform, mapisoform) for m in tmp3]
                            mutations.extend(tmp3)
                            if do_cancer_type:
                                out_metadata['cancer_type'].extend([[cancer_type]]*len(tmp3))
                            if do_cancer_study:
                                out_metadata['cancer_study'].extend([[cancer_study_id]]*len(tmp3))
                            if do_genomic_coordinates:
                                out_metadata['genomic_coordinates'].extend([[None]]*len(tmp3))
                            if do_genomic_mutations:
                                out_metadata['genomic_mutations'].extend([[None]]*len(tmp3))
                else:
                    self.log.debug("no mutations found in this case set")

                self.log.debug("doing case set %s for genetic profiles %s" % (case_set, ",".join(cancer_study_mutation_ids)))
                if len(cancer_study_mutation_ext_ids) > 0:
                    try:
                        response = rq.post(self._requests_url, {'cmd':'getMutationData',
                                                                'case_set_id':case_set,
                                                                'genetic_profile_id':",".join(cancer_study_mutation_ext_ids),
                                                                'gene_list':gene_id}).text
                    except:
                        self.log.error("failed to fetch profile data for these case set and genetic profiles")
                        return None

                    for line in response.strip().split("\n"):
                        if not line.startswith("#") and not line.startswith("entrez_gene_id") and line:
                            tmp=line.strip().split("\t")

                            try:
                                if tmp[5] != "Missense_Mutation":
                                    continue
                            except:
                                pass

                            if mut_prog.match(tmp[7]):
                                mutations.append(self._convert_isoform(tmp[7], mainisoform, cancer_study_isoform, mapisoform))

                            gd = ['hg19', tmp[12], tmp[13], tmp[14], tmp[15]]
                            
                            if tmp[13] != tmp[14]:
                                self.log.warning("mutation corresponds to multiple genomic mutations, genomic mutation won't be annotated")
                                gm = None
                            else:
                                gm = [gd[0], gd[1], gd[2], gd[4], tmp[16]]

                            if do_cancer_type:
                                out_metadata['cancer_type'].append([cancer_type])
                            if do_cancer_study:
                                out_metadata['cancer_study'].append([cancer_study_id])
                            if do_genomic_coordinates:
                                out_metadata['genomic_coordinates'].append(gd)
                            if do_genomic_mutations:
                                out_metadata['genomic_mutations'].append(gm)


        return mutations, out_metadata

    def _convert_isoform(self, mutation, main_isoform, alternate_isoform, mapping):
        if main_isoform == alternate_isoform:
            return mutation
        self.log.debug("converting isoform from %d to %d (%s)" % (main_isoform, mainisoform, cancer_study_id, mutation))
        converted_mutations = []
        resn = int(mutation[1:-1])
        converted = "%s%d%s" % (mutation[0], 
                                mapisoform[cancer_study_isoform][resn],
                                mutation[-1])
        return converted

class COSMIC(DynamicSource, object):
    @logger_init
    def __init__(self, database_files=None, cancer_type=None):
        description = "COSMIC Database"
        super(COSMIC, self).__init__(name='COSMIC', version='v85', description=description)

        if database_files is None:
            database_dir   = '/data/databases/cosmic-v85'
            databases      = [  'CosmicMutantExport',
                                'CosmicCompleteTargetedScreensMutantExport',
                                'CosmicGenomeScreensMutantExport' ]
            self._database_files = [os.path.join(database_dir, i)+'.tsv' for i in databases]
            self._mut_regexp = 'p\.[A-Z][0-9]+[A-Z]$'
            self._mut_prog = re.compile(self._mut_regexp)

    def _get_all_cancer_types(self):
        self.log.info("getting cancer types...")
        cancer_types = set()
        for f in self._database_files:
            with open(f) as fh:
                fh.next()
                for line in fh:
                    cancer_types.add(line.strip().split("\t")[11])
        return sorted(list(cancer_types))

    def _parse_db_file(self, gene_id, cancer_types=[], metadata=[], filter_snps=True):

        mutations = []

        out_metadata = dict(zip(metadata, [list() for i in range(len(metadata))]))

        do_cancer_type = False
        if 'cancer_type' in metadata:
            do_cancer_type = True
        do_genomic_coordinates = False
        if 'genomic_coordinates' in metadata:
            do_genomic_coordinates = True
        do_genomic_mutations = False
        if 'genomic_mutations' in metadata:
            do_genomic_mutations = True


        for f in self._database_files:
            with open(f) as fh:
                fh.next()
                for line in fh:
                    tmp = line.strip().split("\t")
                    if tmp[0] != gene_id or tmp[11] not in cancer_types:
                        continue

                    if filter_snps and tmp[25] == 'y':
                        self.log.info("mutation %s was filtered out as it is classified as SNP" % tmp[18])
                        continue

                    mut_str = tmp[18]
                    if self._mut_prog.match(mut_str):
                        mutations.append(mut_str)
                        if do_cancer_type:
                            out_metadata['cancer_type'].append([tmp[11]])
                        if do_genomic_coordinates or do_genomic_mutations:
                            gd = []
                            if tmp[22] == '38':
                                gd.append('hg38') # genome version
                            elif tmp[22] == '37':
                                gd.append('hg19')
                            else:
                                out_metadata['genomic_coordinates'].append(None)
                                continue
                            tmp2 = tmp[23].split(":")
                            gd.append(tmp2[0]) # chr
                            gd.extend(tmp2[1].split("-")) # [start, end]
                            gd.append(tmp[17][-3]) # ref

                        if do_genomic_coordinates:
                            out_metadata['genomic_coordinates'].append(gd)

                        if do_genomic_mutations:
                            if gd[2] != gd[3]:
                                self.log.warning("mutation corresponds to multiple genomic mutations, genomic mutation won't be annotated")
                                gm = None
                            else:
                                gm = [gd[0], gd[1], gd[2], gd[4], tmp[17][-1]]
                            out_metadata['genomic_mutations'].append(gm)


        return mutations, out_metadata

    def add_mutations(self, sequence, cancer_types=None, use_alias=None, filter_snps=True, metadata=[]):

        if cancer_types is None:
            cancer_types = self._get_all_cancer_types()
            self.log.info("no cancer type specified; will use all of them")
        if use_alias is not None:
            gene_id = sequence.aliases[use_alias]
            self.log.info("using alias %s as gene name" % sequence.aliases[use_alias])
        else:
            gene_id = sequence.gene_id
        raw_mutations, out_metadata = self._parse_db_file(gene_id, cancer_types=cancer_types, metadata=metadata)
        mutations = map(lambda x: x[2:], raw_mutations)
        unique_mutations = list(set(mutations))
        self.log.info("unique mutations found in %s: %s" % (self.name, ", ".join(sorted(unique_mutations, key=lambda x: int(x[1:-1])))))

        for m in unique_mutations:
            wt = m[0]
            mut = m[-1]
            num = int(m[1:-1])

            if wt == mut:
                self.log.info("synonymous mutation %s discarded" % m)
                continue

            try:
                site_seq_idx = sequence.seq2index(num)
            except:
                self.log.warning("mutation %s is outside the protein sequence; it will be skipped")
                continue

            position = sequence.positions[site_seq_idx]

            if position.wt_residue_type != wt:
                self.log.warning("for mutation %s, residue %d is %s in wild-type sequence; it will be skipped" %(m, wt, position.wt_residue_type))
                continue

            mutation_indices = [i for i, x in enumerate(mutations) if x == m]

            mutation_obj = Mutation(sequence.positions[site_seq_idx],
                                    mut,
                                    [self])

            for md in metadata:
                mutation_obj.metadata[md] = []
                for mi in mutation_indices:
                    if out_metadata[md][mi] is not None:
                        tmp_md = [self] + out_metadata[md][mi]
                        this_md = metadata_classes[md](*tmp_md)
                        mutation_obj.metadata[md].append(this_md)
            position.add_mutation(mutation_obj)

class PhosphoSite(DynamicSource, object):
    @logger_init
    def __init__(self, database_file=None):
        description = "PhosphoSite Database"
        super(PhosphoSite, self).__init__(name='PhosphoSite', version='1.0', description=description)
        if database_file is None:
            self._ptm_types = ['acetylation', 'methylation', 'O-GalNAc', 'O-GlcNAc', 'phosphorylation', 'sumoylation', 'ubiquitination']
            self._ptm_suffixes = ['ac', 'm[0-9]', 'ga', 'gl', 'p', 'sm', 'ub']
            self._ptm_suffix_offsets = [-3, -3, -3, -3, -2, -3, -3]
            self._database_dir = "/data/databases/phosphosite/"
            self._database_files = dict([(i, "%s/%s_site_dataset"%(self._database_dir,i)) for i in self._ptm_types])

    def _parse_db_file(self, gene_id):

        sites = dict(zip(self._ptm_types, [list() for i in range(len(self._ptm_types))]))

        for ptm_idx,ptm in enumerate(self._ptm_types):
            p_regexp = '[A-Z][0-9]+-%s' % self._ptm_suffixes[ptm_idx]
            p_prog = re.compile(p_regexp)
            p_sites = []

            with open(self._database_files[ptm]) as fh:
                fh.next()
                fh.next()
                fh.next()
                for line in fh:
                    tmp = line.strip().split("\t")
                    if str.upper(tmp[0]) != gene_id or tmp[6] != 'human':
                        continue
                    p_str = tmp[4]
                    if p_prog.match(p_str):
                        p_sites.append(p_str[:self._ptm_suffix_offsets[ptm_idx]])
            sites[ptm] = p_sites
        return sites

    def add_position_properties(self, sequence):
        sites = self._parse_db_file(sequence.gene_id)
        for ptm_idx,ptm in enumerate(self._ptm_types):

            p_sites = sites[ptm]
            unique_p_sites = list(set(p_sites))
            for m in unique_p_sites:
                wt = m[0]
                site = int(m[1:])

                try:
                    site_seq_idx = sequence.seq2index(num)
                except:
                    self.log.warning("mutation %s is outside the protein sequence; it will be skipped")
                    continue

                position = sequence.positions[site_seq_idx]
                if position.wt_residue_type != wt:
                    self.log.warning("for mutation %s, residue %d is %s in wild-type sequence; it will be skipped" %(m, wt, position.wt_residue_type))
                    continue

                already_annotated = False
                for prop in position.properties:
                    if isinstance(prop, position_properties_classes[ptm]):
                        prop.sources.append(self)
                        self.log.info("site %s already annotated as %s; source will be added" % (m, position_properties_classes[ptm].name))
                        already_annotated = True
                
                if not already_annotated:
                    property_obj = position_properties_classes[ptm](  sources=[self],
                                                        position=sequence.positions[site_seq_idx]
                                                        )
                    position.add_property(property_obj)
                    self.log.info("adding %s to site %s" % (m, property_obj.name))

class MyVariant(DynamicSource, object):
    @logger_init
    def __init__(self):
        description = "MyVariant.Info Database"
        super(MyVariant, self).__init__(name='MyVariant', version='1,0', description=description)

        self._mv = myvariant.MyVariantInfo()
        self._lo = pyliftover.LiftOver('hg38', 'hg19')

    def add_metadata(self, sequence, md_type='revel'):
        if md_type == 'revel':
            add_this_metadata = self._get_revel

        for pos in sequence.positions:
            for mut in pos.mutations:
                add_this_metadata(mut)

    def _get_revel(self, mutation):
        if 'genomic_coordinates' not in mutation.metadata.keys():
            self.log.warning("no genomic coordiantes for mutation %s; it will be skipped" % mutation)
            return False

        mutation.metadata['revel_score'] = list()

        gcs = list(set(mutation.metadata['genomic_coordinates']))

        for gc in gcs:
            if gc.coord_start != gc.coord_end:
                self.log.warning("mutation %s has more than one nucleotide change; it will be skipped" % mutation)
                continue
            if gc.genome_version == 'hg38':
                self.log.info("mutation %s has genomic data in hg38 assembly - will be converted to hg19" % mutation)
                converted_coords = self._lo.convert_coordinate('chr%s' % gc.chr, int(gc.coord_start))
                assert len(converted_coords) == 1
                converted_coords = (converted_coords[0][0], converted_coords[0][1])
            elif gc.genome_version == 'hg19':
                converted_coords = ('chr%s' % gc.chr, int(gc.coord_start))
            else:
                self.log.error("genomic coordinates are not expressed either in hg38 or hg19 for mutation %s; it will be skipped" % mutation)
                return False

            query_str = 'chr%s:%d' % converted_coords
            query = self._mv.query(query_str)
            hits = query['hits']
            if len(hits) < 1:
                self.log.warning("no DBnsfp hits for mutation %s; it will be skipped" % mutation)
                continue
            else:
                self.log.info("%d DBnsfp hits for mutation %s" % (len(hits), mutation))
            revel_scores = []
            for hit in hits:
                try:
                    hit_pos = hit['dbnsfp']['aa']['pos']
                    hit_ref = hit['dbnsfp']['aa']['ref']
                    hit_alt = hit['dbnsfp']['aa']['alt']
                except KeyError:
                    self.log.warning("no residue information was found in DBsnfp hit. It will be skipped")
                    continue

                if not mutation.sequence_position.wt_residue_type == hit['dbnsfp']['aa']['ref']:
                    self.log.warning("reference residue in revel does not correspond; it will be skipped")
                    continue

                if not mutation.sequence_position.sequence_position == hit_pos:
                    try:
                        if not mutation.sequence_position.sequence_position in map(int, hit['dbnsfp']['aa']['pos']):
                            raise TypeError
                    except:
                        self.log.warning("sequence position in revel does not correspond for this hit; it will be skipped")
                        continue

                if mutation.mutated_residue_type != hit_alt:
                    self.log.info("protein mutation in revel does not correspond for this hit; it will be skipped")
                    continue

                revel_score = None
                try:
                    revel_score = hit['dbnsfp']['revel']['score']
                except KeyError:
                    self.log.warning("no revel score found for mutation %s in this hit; it will be skipped" % mutation)
                    continue

                revel_scores.append(revel_score)

            revel_scores = list(set(revel_scores))
            if len(revel_scores) > 1:
                self.log.warning("more than one revel score found; it will be skipped (%s)" % ", ".join(['%.3f' % i for i in revel_scores]))
                return False
            elif len(revel_scores) == 0:
                self.log.warning("no revel found at all for mutation %s" % mutation)
                return False
            else:
                revel_md = DbnsfpRevel(source=self, score=revel_score)

                mutation.metadata['revel_score'].append(revel_md)

        mutation.metadata['revel_score'] = list(set(mutation.metadata['revel_score']))
        self.log.info("found revel scores for mutation %s: %s" % (mutation, ", ".join(["%.3f"%i.score for i in mutation.metadata['revel_score']])))
        return True

class ELMDatabase(DynamicSource, object):
    def __init__(self):
        description = "ELM Database"
        super(MyVariant, self).__init__(name='ELM', version='1.0', description=description)

        self._requests_url = "http://elm.eu.org/elms/"
        self._get_elm_classes()

    def _get_elm_classes(self):
        self._elm_classes = {}

        response = rq.get(self._requests_url+'/elms_index.tsv').text
        tmp = response.split("\n")
        for line in tmp:
            if line.startswith('"ELME'):
                tmp2 = line.strip().split()
                self._elm_classes[tmp2[0]] = tmp2[1:3]

    def _get_annotations(self, gene_name):
        pass

class ELMPredictions(DynamicSource, object):
    @logger_init
    def __init__(self):
        description = "ELM Prediction"
        super(ELMPredictions, self).__init__(name='ELM', version='1.0', description=description)

        self._classes_url  = "http://elm.eu.org/elms/"
        self._requests_url = "http://elm.eu.org/start_search/"
        self._get_elm_classes()

    def _get_elm_classes(self):
        self._elm_classes = {}
        unquote = lambda x: str.strip(x, '"')
        self.log.info("retrieving ELM classes")
        try:
            response = rq.get(self._classes_url+'/elms_index.tsv').text
        except:
            self.log.error("couldn't retrieve ELM classes")
            return
        tmp = response.split("\n")
        for line in tmp:
            if line.startswith('"ELME'):
                tmp2 = line.strip().split("\t")
                tmp2 = map(str, tmp2)
                tmp2 = map(unquote, tmp2)
                self._elm_classes[tmp2[1]] = tmp2[2:]

    def _get_prediction(self, gene_name):
        self.log.info("retrieving prediction for %s" % gene_name )
        try:
            req_url = os.path.join(self._requests_url, gene_name) + ".tsv"
        except:
            self.log.error("couldn't fetch ELM predictions")
            return None

        out = []
        response = rq.get(req_url).text
        tmp = response.split("\n")
        for line in tmp:
            tmp2 = line.strip().split()
            if line.startswith("#") or line.startswith('elm_identifier') or not line:
                continue
            assert len(tmp2) == 10
            tmp2[0] = str(tmp2[0])
            tmp2[1] = int(tmp2[1])
            tmp2[2] = int(tmp2[2])
            for i in range(len(tmp2[3:])):
                tmp2[i+3] = tmp2[i+3] == 'True'
            out.append(tmp2)
        return out

    def add_sequence_properties(self, sequence, exclude_elm_classes='{100}', use_alias=None):
        self.log.info("adding ELM predictions to sequence ...")
        if use_alias is None:
            data = self._get_prediction(sequence.gene_id)
        else:
            self.log.info("will use alias %s as gene name" % sequence.aliases[use_alias])
            data = self._get_prediction(sequence.aliases[use_alias])

        for d in data:
            if d[5]:
                self.log.info("%s was filtered out by ELM")
                continue

            if re.match(exclude_elm_classes, d[0]):
                self.log.info("%s was filtered out as requested" % d[0])
                continue

            this_positions = []
            for p in range(d[1],d[2]+1):
                this_positions.append(sequence.positions[sequence.seq2index(p)])

            property_obj = sequence_properties_classes['linear_motif']  (sources=[self],
                                                                         positions=this_positions,
                                                                         lmtype=self._elm_classes[d[0]][0])

            property_obj.metadata['function'] = [self._elm_classes[d[0]][0]]
            property_obj.metadata['ref']      = self.description
            sequence.add_property(property_obj)

class ManualAnnotation(StaticSource):
    @logger_init
    def __init__(self, datafile, csvoptions=None):
        description="Annotations from %s" % datafile
        super(ManualAnnotation, self).__init__(name='Manual annotations', version='', description=description)

        self._datafile = datafile
        fh = open(self._datafile, 'rb')
        if csvoptions is None:
            csvoptions = {}
        try:
            self._csv = csv.reader(fh, **csvoptions)
        except:
            log.error("")
        self._table = None

    def _parse_datafile(self):
        self.log.info("Parsing annotations from %s" % self._datafile)
        out = [] # type seq type function ref
        ptm_keywords = ['cleavage', 'phosphorylation', 'ubiquitination', 'acetylation', 'sumoylation', 's-nitrosylation', 'methylation']
        for idx,row in enumerate(self._csv):
            row_number = idx+1
            if len(row) < 5:
                self.log.warning("row %d has less than 5 columns and will be skipped. You need at least 5 columns to be present." % row_number)
                continue
            if   row[2] in ptm_keywords:
                this_row = [row[2], tuple([int(row[1])]), row[0], row[3], row[4]]
                self.log.info("added %s from row %d" % (row[2], row_number))
            elif row[2] == 'linear_motif':
                this_row = [row[2]]
                if '-' in row[1]:
                    tmp = map(int, row[1].split("-"))
                    this_row.append(tuple(range(tmp[0], tmp[1]+1)))
                else:
                    this_row.append( (int(row[1]) ))
                self.log.info("added %s from row %d" % (row[2], row_number))
                this_row.extend([row[0], row[3], row[4]])
            elif row[2] == 'mutation':
                this_row = [ row[2], tuple([int(row[1][1:-1])]), row[1], row[0], row[3], row[4] ]
                self.log.info("added %s from row %d" % (row[2], row_number))
            else:
                self.log.warning("no recognized keyword was found in column 3 and row %d will be skipped." % row_number)
                continue

            out.append(tuple(this_row))
        self._table = out


    def add_mutations(self, sequence):
        if self._table is None:
            self._parse_datafile()
        else:
            self.log.info("datafile information has already been parsed - won't be re-read")

        mutations = [d[2] for d in self._table if d[0] == 'mutation']
        unique_mutations = list(set(mutations))

        self.log.info("unique mutations found in datafile: %s" % (", ".join(sorted(unique_mutations, key=lambda x: int(x[1:-1])))))

        for m in unique_mutations:
            wt = m[0]
            mut = m[-1]
            num = int(m[1:-1])

            if wt == mut:
                self.log.info("synonymous mutation %s discarded" % m)
                continue

            try:
                site_seq_idx = sequence.seq2index(num)
            except:
                self.log.warning("mutation %s is outside the protein sequence; it will be skipped")
                continue

            position = sequence.positions[site_seq_idx]

            if position.wt_residue_type != wt:
                self.log.warning("for mutation %s, residue %d is %s in wild-type sequence; it will be skipped" %(m, wt, position.wt_residue_type))
                continue

            mutation_indices = [i for i, x in enumerate(mutations) if x == m]

            mutation_obj = Mutation(sequence.positions[site_seq_idx],
                                    mut,
                                    [self])

            position.add_mutation(mutation_obj)

    def add_position_properties(self, sequence):
        if self._table is None:
            self._parse_datafile()
        else:
            self.log.info("datafile information has already been parsed - won't be re-read")

        for d in self._table:
            this_position = sequence.positions[sequence.seq2index(d[1][0])]
            
            try:
                property_obj = position_properties_classes[d[0]](sources=[self],
                                                                 position=this_position)
                assert len(d[1]) == 1
            except KeyError:
                self.log.warning("property %s is not a valid position property; row will be skipped" % d[0])
                continue

            property_obj.metadata['function'] = [d[3]]
            property_obj.metadata['ref']      = [d[4]]
            this_position.add_property(property_obj)

    def add_sequence_properties(self, sequence):
        if self._table is None:
            self._parse_datafile()
        else:
            self.log.info("datafile information has already been parsed - won't be re-read")

        for d in self._table:
            this_positions = []
            for p in d[1]:
                this_positions.append(sequence.positions[sequence.seq2index(p)])

            try:
                property_obj = sequence_properties_classes[d[0]](sources=[self],
                                                             positions=this_positions,
                                                             lmtype=d[2])
            except KeyError:
                self.log.warning("property %s is not a valid sequence property; row will be skipped" % d[0])
                continue

            property_obj.metadata['function'] = [d[3]]
            property_obj.metadata['ref']      = [d[4]]
            sequence.add_property(property_obj)
