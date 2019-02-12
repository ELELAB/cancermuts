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
import pandas as pd
from .core import Sequence, Mutation
from .properties import *
from .metadata import *
from .log import *
import json
import re
import sys
import os
import csv
import myvariant
import pyliftover
from future.utils import iteritems

import sys
if sys.version_info[0] >= 3:
    unicode = str


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

    default_strand = '+'
    
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
                self.log.warning("for mutation %s, residue %s is %s in wild-type sequence; it will be skipped" %(m, wt, position.wt_residue_type))
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
            data.append(tuple(map(unicode, line.strip().split("\t"))))
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
        for cancer_study_id in list(zip(*self._cache_cancer_studies))[0]:
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
        for cancer_study_id in list(zip(*self._cache_cancer_studies))[0]:
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

        out_metadata = dict(list(zip(metadata, [list() for i in range(len(metadata))])))

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
                    self.log.error("cancer type for %s not found - cancer study will be used instead" % cancer_study_suffix)
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
                            tmp2 = [x for x in tmp if x != 'NaN']
                            tmp3 = [x for x in tmp2 if mut_prog.match(x)]
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

                            if do_cancer_type:
                                out_metadata['cancer_type'].append([cancer_type])
                            if do_cancer_study:
                                out_metadata['cancer_study'].append([cancer_study_id])
                            if do_genomic_coordinates:
                                gd = ['hg19', tmp[12], tmp[13], tmp[14], tmp[15]]
                                out_metadata['genomic_coordinates'].append(gd)
                            if do_genomic_mutations:
                                if tmp[13] != tmp[14]:
                                    self.log.warning("mutation corresponds to multiple genomic mutations, genomic mutation won't be annotated")
                                    gm = None
                                else:
                                    gm = ['hg19', tmp[12], self.default_strand, gd[2], gd[4], tmp[16]]
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
        super(COSMIC, self).__init__(name='COSMIC', version='v87', description=description)

        self._mut_regexp = 'p\.[A-Z][0-9]+[A-Z]$'
        self._mut_prog = re.compile(self._mut_regexp)
        dataframes = []

        if database_files is None:
            database_dir   = '/data/databases/cosmic'
            databases      = [  'CosmicMutantExport',
                                'CosmicCompleteTargetedScreensMutantExport',
                                'CosmicGenomeScreensMutantExport' ]
            self._database_files = [os.path.join(database_dir, i)+'.tsv' for i in databases]
        else:
            self._database_files = database_files

        for f in self._database_files:
            try:
                self.log.info("Parsing database file %s ..." % f)
                dataframes.append( pd.read_csv(f, sep='\t', dtype='str', na_values='NS') )
            except:
                self.log.error("Couldn't parse database file.")
                self._df = None
                return

        self._df = pd.concat(dataframes, ignore_index=True, sort=False)

    def _parse_db_file(self, gene_id, cancer_types=None, metadata=[], filter_snps=True):

        site_kwd = ['Primary_site', 'Site_subtype_1', 'Site_subtype_2', 'Site_subtype_3']
        histology_kwd = ['Primary_histology', 'Histology_subtype_1', 'Histology_subtype_2', 'Histology_subtype_3']

        mutations = []

        out_metadata = dict(list(zip(metadata, [list() for i in range(len(metadata))])))

        do_cancer_type = False
        if 'cancer_type' in metadata:
            do_cancer_type = True
        do_genomic_coordinates = False
        if 'genomic_coordinates' in metadata:
            do_genomic_coordinates = True
        do_genomic_mutations = False
        if 'genomic_mutations' in metadata:
            do_genomic_mutations = True
        do_site = False
        if 'cancer_site' in metadata:
            out_metadata['cancer_site'] = []
            do_site = True
        do_histology = False
        if 'cancer_histology' in metadata:
            out_metadata['cancer_histology'] = []
            do_histology = True

        df = self._df[ self._df['Gene name'] == gene_id ]

        if filter_snps:
            df = df[ df['SNP'] != 'y' ]

        if cancer_types is not None:
            df = df[ df['Primary histology'] in cancer_types ]

        df = df[ df['Mutation AA'].notna() ]

        df = df[ df.apply(lambda x: bool(self._mut_prog.match(x['Mutation AA'])), axis=1) ]

        df.columns = df.columns.str.replace(r'\s+', '_')

        for r in df.itertuples():
            mutations.append(r.Mutation_AA)

            if do_cancer_type:
                out_metadata['cancer_type'].append([r.Primary_histology])

            if do_genomic_coordinates or do_genomic_mutations:
                gd = []
                grch = r.GRCh
                if grch == '38':
                    gd.append('hg38') # genome version
                elif grch == '37':
                    gd.append('hg19')
                else:
                    out_metadata['genomic_coordinates'].append(None)
                    self.log.warning("Genome assembly not specified for mutation; genomic coordinates won't be annotated")
                    continue

                tmp2 = r.Mutation_genome_position.split(":")
                gd.append(tmp2[0]) # chr
                gd.extend(tmp2[1].split("-")) # [start, end]
                gd.append(r.Mutation_CDS[-3]) # ref

            if do_genomic_coordinates:
                out_metadata['genomic_coordinates'].append(gd)

            if do_genomic_mutations:
                if gd[2] != gd[3]:
                    self.log.warning("mutation corresponds to multiple genomic mutations, genomic mutation won't be annotated")
                    gm = None
                else:
                    gm = [gd[0], gd[1], r.Mutation_strand, gd[2], gd[4], r.Mutation_CDS[-1]]
                out_metadata['genomic_mutations'].append(gm)

            if do_site:
                out_metadata['cancer_site'].append([r.__getattribute__(a) for a in site_kwd])

            if do_histology:
                out_metadata['cancer_histology'].append([r.__getattribute__(a) for a in histology_kwd])

        return mutations, out_metadata

    def add_mutations(self, sequence, cancer_types=None, use_alias=None, filter_snps=True, metadata=[]):

        if cancer_types is None:
            self.log.info("no cancer type specified; will use all of them")
        if use_alias is not None:
            gene_id = sequence.aliases[use_alias]
            self.log.info("using alias %s as gene name" % sequence.aliases[use_alias])
        else:
            gene_id = sequence.gene_id

        raw_mutations, out_metadata = self._parse_db_file(gene_id, cancer_types=cancer_types, metadata=metadata)
        mutations = [x[2:] for x in raw_mutations]
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
                self.log.warning("for mutation %s, residue %s is %s in wild-type sequence; it will be skipped" %(m, wt, position.wt_residue_type))
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
    def __init__(self, database_files=None):
        description = "PhosphoSite Database"
        super(PhosphoSite, self).__init__(name='PhosphoSite', version='1.0', description=description)

        self._ptm_types = ['acetylation', 'methylation', 'O-GalNAc', 'O-GlcNAc', 'phosphorylation', 'sumoylation', 'ubiquitination']
        self._ptm_types_to_classes = {  'acetylation'     : 'ptm_acetylation', 
                                        'methylation'     : 'ptm_methylation', 
                                        'O-GalNAc'        : 'ptm_ogalnac', 
                                        'O-GlcNAc'        : 'ptm_oglcnac', 
                                        'phosphorylation' : 'ptm_phosphorylation', 
                                        'sumoylation'     : 'ptm_sumoylation',
                                        'ubiquitination'  : 'ptm_ubiquitination' }
        self._ptm_suffixes = ['ac', 'm[0-9]', 'ga', 'gl', 'p', 'sm', 'ub']
        self._ptm_suffix_offsets = [-3, -3, -3, -3, -2, -3, -3]

        if database_files is None:
            self._database_dir = "/data/databases/phosphosite/"
            self._database_files = dict([(i, "%s/%s_site_dataset"%(self._database_dir,i)) for i in self._ptm_types])

        else:
            self._database_files = database_files

        self._dataframes = {}
        for k,f in iteritems(self._database_files):
            try:
                self._dataframes[k] = pd.read_csv(f, skiprows=3, sep='\t')
                self._dataframes[k] = self._dataframes[k][ self._dataframes[k]['ORGANISM'] == 'human']
            except:
                self.log.error("couldn't read database file %s" % f)



    def _parse_db_file(self, gene_id):

        sites = dict(list(zip(self._ptm_types, [list() for i in range(len(self._ptm_types))])))

        for ptm_idx,ptm in enumerate(self._ptm_types):
            p_regexp = '[A-Z][0-9]+-%s' % self._ptm_suffixes[ptm_idx]
            p_prog = re.compile(p_regexp)
            p_sites = []

            df = self._dataframes[ptm]
            df = df[ df['GENE'].str.upper() == gene_id.upper() ]
            df = df[ df.apply(lambda x: bool(p_prog.match(x['MOD_RSD'])), axis=1) ]

            sites[ptm] = [ x[:self._ptm_suffix_offsets[ptm_idx]] for x in df['MOD_RSD'].values  ]

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
                    site_seq_idx = sequence.seq2index(site)
                except:
                    self.log.warning("PTM site %s is outside the protein sequence; it will be skipped" % m)
                    continue

                position = sequence.positions[site_seq_idx]
                if position.wt_residue_type != wt:
                    self.log.warning("for PTM %s, residue %d is %s in wild-type sequence; it will be skipped" %(m, wt, position.wt_residue_type))
                    continue

                already_annotated = False
                for prop in position.properties:
                    if isinstance(prop, position_properties_classes[ptm]):
                        prop.sources.append(self)
                        self.log.info("site %s already annotated as %s; source will be added" % (m, position_properties_classes[ptm].name))
                        already_annotated = True
                
                if not already_annotated:
                    property_obj = position_properties_classes[self._ptm_types_to_classes[ptm]](  sources=[self],
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

        self._supported_metadata = {'revel_score' : self._get_revel}


    def add_metadata(self, sequence, md_type=['revel_score']):
        if type(md_type) is str:
            md_types = [md_type]
        else:
            md_types = md_type

        metadata_functions = []

        for md_type in md_types:
            try:
                metadata_functions.append(self._supported_metadata[md_type])
            except KeyError:
                self.log.warning("MyVariant doesn't support metadata type %s" % md_type)

        for pos in sequence.positions:
            for mut in pos.mutations:
                for add_this_metadata in metadata_functions:
                    add_this_metadata(mut)


    def _get_revel(self, mutation):

        revel_score = None

        mutation.metadata['revel_score'] = list()

        if 'genomic_coordinates' not in mutation.metadata and 'genomic_mutations' not in mutation.metadata:
            self.log.warning("no genomic coordinates or genomic mutations data for mutation %s; it will be skipped" % mutation)
            return False

        if 'genomic_coordinates' in mutation.metadata:
            gcs = mutation.metadata['genomic_coordinates']
        else:
            self.log.warning("no genomic coordinates available for revel score, mutation %s" % mutation)
            gcs = [None] * len(mutation.metadata['genomic_mutations'])

        if 'genomic_mutations' in mutation.metadata:
            gms = mutation.metadata['genomic_mutations']
        else:
            self.log.warning("no genomic mutations available for revel score, mutation %s" % mutation)
            gms = [None] * len(mutation.metadata['genomic_coordinates'])

        gcs_gms = list(zip(gcs, gms))

        for gc,gm in gcs_gms:
            if gm is not None:
                self.log.info("genomic mutation will be used to retrieve revel score for mutation %s" % mutation)
                revel_score = self._get_revel_from_gm(mutation, gm)
            if gc is not None and revel_score is None:
                self.log.info("genomic coordinates will be used to retrieve revel score for mutation %s" % mutation)
                revel_score = self._get_revel_from_gc(mutation, gc)
            if revel_score is None:
                self.log.warning("mutations %s has no genomic coordinates or genomic mutation available; it will be skipped" % mutation)
                return

            mutation.metadata['revel_score'].append(revel_score)

    def _validate_revel_hit(self, mutation, hit):

        try:
            aa = hit['dbnsfp']['aa']
        except:
            self.log.error("no residue information was found; it will be skipped")
            return False

        if isinstance(aa, list):
            if len(aa) == 1:
                aa = aa[0]
            else:
                self.log.error("more than one aminoacid specified; it will be skipped")
                return False
        try:
            hit_pos = aa['pos']
            hit_ref = aa['ref']
            hit_alt = aa['alt']
        except:
            self.log.warning("no residue information was found in the variant information; it will be skipped")
            return False
    
        if not mutation.sequence_position.wt_residue_type == aa['ref']:
            self.log.warning("reference residue in revel does not correspond; it will be skipped")
            return False

        if not mutation.sequence_position.sequence_position == hit_pos:
            try:
                if not mutation.sequence_position.sequence_position in [ int(a) for a in aa['pos'] ]:
                    raise TypeError
            except:
                self.log.warning("sequence position in revel does not correspond for this hit; it will be skipped")
                return False

        if mutation.mutated_residue_type != hit_alt:
            self.log.info("protein mutation in revel does not correspond for this hit; it will be skipped")
            return False

        return True

    def _convert_hg38_to_hg19(self, gc):
        if gc.genome_version == 'hg38':
            self.log.info("%s has genomic data in hg38 assembly - will be converted to hg19" % gc)
            converted_coords = self._lo.convert_coordinate('chr%s' % gc.chr, int(gc.get_coord()))
            assert len(converted_coords) == 1
            converted_coords = (converted_coords[0][0], converted_coords[0][1])
        elif gc.genome_version == 'hg19':
            converted_coords = ('chr%s' % gc.chr, int(gc.get_coord()))
        else:
            self.log.error("genomic coordinates are not expressed either in hg38 or hg19 for %s; it will be skipped" % gc)
            converted_coords = None
        return converted_coords

    def _get_revel_from_gm(self, mutation, gm):

        converted_coords = self._convert_hg38_to_hg19(gm)
        if converted_coords is None:
            return None

        query_str = '%s:g.%d%s>%s' % (converted_coords[0], converted_coords[1], gm.get_sense_wt(), gm.get_sense_mut())

        hit = self._mv.getvariant(query_str)

        if hit is None:
            self.log.error("variant %s was not found in MyVariant!" % query_str)
            return None

        if not self._validate_revel_hit(mutation, hit):
            return None
        
        try:
            revel_score = hit['dbnsfp']['revel']['score']
        except KeyError:
            self.log.warning("no revel score found for mutation %s in this hit; it will be skipped" % mutation)
            return None

        return DbnsfpRevel(source=self, score=revel_score)

    def _get_revel_from_gc(self, mutation, gc):

        found_scores = []

        if gc is None:
            return None

        if gc.coord_start != gc.coord_end:
            self.log.warning("mutation %s has more than one nucleotide change; it will be skipped" % mutation)
            return None

        converted_coords = self._convert_hg38_to_hg19(gc)
        if converted_coords is None:
            return None

        query_str = 'chr%s:%d' % converted_coords
        query = self._mv.query(query_str)
        hits = query['hits']
        if len(hits) < 1:
            self.log.warning("no DBnsfp hits for mutation %s; it will be skipped" % mutation)
            return None
        else:
            self.log.info("%d DBnsfp hits for mutation %s" % (len(hits), mutation))
        revel_scores = []
 
        for hit in hits:

            revel_score = None

            if not self._validate_revel_hit(mutation, hit):
                continue

            try:
                revel_score = hit['dbnsfp']['revel']['score']
            except KeyError:
                self.log.warning("no revel score found for mutation %s in this hit; it will be skipped" % mutation)
                return None

            revel_scores.append(revel_score)

        revel_scores = list(set(revel_scores))
        if len(revel_scores) > 1:
            self.log.warning("more than one revel score found; it will be skipped (%s)" % ", ".join(['%.3f' % i for i in revel_scores]))
            return None
        elif len(revel_scores) == 0:
            self.log.warning("no revel scores found using genomic coordinates for mutation %s" % mutation)
            return None
        else:
            revel_md = DbnsfpRevel(source=self, score=revel_score)
            return revel_md

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
                tmp2 = [ unicode(t) for t in tmp2 ]  
                tmp2 = [ unicode.strip(t, '"') for t in tmp2 ]
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

class ExAC(DynamicSource, object):
    
    description = "ExAC"

    @logger_init
    def __init__(self):

        super(ExAC, self).__init__(name='ExAC', version='0.1', description=self.description)

        self._requests_url = 'http://exac.broadinstitute.org/variant'
        self._data_cache = {}
        self._supported_metadata = {'exac_allele_frequency' : self._get_allele_freq,
                                    'exac_af_filter'        : self._get_af_filter}

    def add_metadata(self, sequence, md_type=['exac_allele_frequency','exac_af_filter']):
        
        self.log.debug("Adding metadata: " + ', '.join(md_type) )

        if type(md_type) is str:
            md_types = [md_type]
        else:
            md_types = md_type

        metadata_functions = []

        for md_type in md_types:
            try:
                metadata_functions.append(self._supported_metadata[md_type])
            except KeyError:
                self.log.warning("ExAC doesn't support metadata type %s" % md_type)

        self.log.debug("collected metadata functions: %s" % ', '.join([i for i in metadata_functions.__repr__()]))

        for pos in sequence.positions:
            for mut in pos.mutations:
                for i,add_this_metadata in enumerate(metadata_functions):
                    mut.metadata[md_types[i]] = []
    
                if not 'genomic_mutations' in mut.metadata:
                    self.log.warning("no genomic mutation data available for Exac SNP, mutation %s. It will be skipped" % mut)
                    continue
                elif mut.metadata['genomic_mutations'] is None or len(mut.metadata['genomic_mutations']) == 0:
                    self.log.warning("no genomic mutation data available for Exac SNP, mutation %s. It will be skipped" % mut)
                    continue

                for i,add_this_metadata in enumerate(metadata_functions):
                    self.log.debug("adding metadata %s to %s" % (md_types[i], mut))
                    add_this_metadata(mut)

    def _get_allele_freq(self, mutation):
        self.log.debug("getting allele frequency")
        self._get_metadata(mutation, md_type='exac_allele_frequency')

    def _get_af_filter(self, mutation):
        self.log.debug("getting af filter")
        self._get_metadata(mutation, md_type='exac_af_filter')


    def _get_metadata(self, mutation, md_type):

        exac_key = {    'exac_allele_frequency' : 'allele_freq', 
                        'exac_af_filter'        : 'af_filter'       }

        mutation.metadata[md_type] = list()

        gcs = mutation.metadata['genomic_mutations']

        for gc in gcs:
            try:
                exac = self._data_cache[gc]
            except KeyError:
                self._data_cache[gc] = self._get_exac_data(gc)

            if self._data_cache[gc] is None:
                #self.log.warning("data for variant %s couldn't be retrieved" % gc)
                continue
            
            exac = self._data_cache[gc]

            mutation.metadata[md_type].append(metadata_classes[md_type](self, exac[exac_key[md_type]]))

    def _parse_exac_page(self, text):
        for line in text.split('\n'):
            if line.strip().startswith('window.variant = '):
                data = json.loads(line.split(' = ')[1][:-1])
                return data

    def _get_exac_data(self, variant):
        if type(variant) is GenomicMutation:
            v_str = variant.as_hg19().get_value_str(fmt='exac')
        else:
            v_str = variant

        self.log.debug("variant string is %s" % v_str)  

        try:
            self.log.info("fetching %s" % '/'.join( [ self._requests_url, v_str ] ) )
            request = rq.get('/'.join( [ self._requests_url, v_str ] ))
        except:
            self.log.error("Couldn't perform request for EXaC variant")
            return None

        if request.status_code == 404:
            self.log.warning("the specified variant %s wasn't found on ExAC" % v_str)
            return None
        elif request.status_code == 200:
            try:
                data = self._parse_exac_page(request.text)
            except:
                self.log.warning("Couldn't parse downloaded file for variant %s" % v_str)
                return None
        else:
            self.log.error("The request for EXaC variant returneds an unexpected status code")
            return None

        return data

class MobiDB(DynamicSource):

    description = "MobiDB"

    @logger_init
    def __init__(self):

        super(MobiDB, self).__init__(name='MobiDB', version='0.1', description=self.description)

        self._requests_url = 'http://mobidb.bio.unipd.it/ws/'
        self._data_cache = {}
        self._supported_properties = { 'mobidb_disorder_propensity' : self._get_mobidb_disorder_predictions }
                                    
    def add_position_properties(self, sequence, prop=['mobidb_disorder_propensity'], use_alias=None):

        if type(prop) is str:
            props = [prop]
        else:
            props = prop

        self.log.debug("Adding property: " + ', '.join(props) )

        prop_functions = []

        for prop in props:
            try:
                prop_functions.append(self._supported_properties[prop])
            except KeyError:
                self.log.warning("MobiDB doesn't support property type %s" % prop)

        self.log.debug("collected property functions: %s" % ', '.join([i for i in prop_functions.__repr__()]))

        for i,add_this_property in enumerate(prop_functions):
            self.log.debug("adding property %s to %s" % (props[i], sequence))
            add_this_property(sequence, use_alias=use_alias)


    def _get_mobidb_disorder_predictions(self, sequence, *args, **kwargs):

        assignments = self._get_mobidb_disorder_predictions_assignments(sequence, *args, **kwargs)

        if len(assignments) != len(sequence.positions): 
            self.log.error("Error in the assignment of predictions")
            return False

        for i,a in enumerate(assignments):
            sequence.positions[i].properties['mobidb_disorder_propensity'] = position_properties_classes['mobidb_disorder_propensity'](sequence.positions[i], [self], a)


    def _get_mobidb_disorder_predictions_assignments(self, sequence, *args, **kwargs):
        if 'use_alias' in kwargs:
            use_alias = kwargs['use_alias']
        else:
            use_alias = None

        data = self._get_mobidb(sequence, use_alias=use_alias)

        try:
            disorder_data = data['mobidb_consensus']['disorder']
        except KeyError:
            self.log.error("No disorder data was found")
            return None

        assignments = [None for p in sequence.positions]

        if 'full' in disorder_data:
            self.log.info("full consensus available")

            if len(data['mobidb_consensus']['disorder']['full']) > 0:
                self.log.warning("More than one prediction found; will use the first")

            regions = data['mobidb_consensus']['disorder']['full'][0]['regions']

            for r in regions:
                for i in range(r[0]-1, r[1]):
                    try:
                        assignments[i] = str(r[2])
                    except IndexError:
                        self.log.error("residue index %s not in sequence!" % i)
                        return None
            if None in assignments:
                self.log.warning("Some positions were not assigned!")

        elif 'predictors' in disorder_data:
            self.log.info("full consensus not available - predictions will be used")

            mobidb_lite = data['mobidb_consensus']['disorder']['predictors'][1]
            if mobidb_lite['method'] != 'mobidb-lite':
                self.log.error("ModiDB-lite consensus prediction is not available")
                return None

            regions = mobidb_lite['regions']
            for r in regions:
                for i in range(r[0]-1,r[1]): # r[1]: -1 because of the 0-offset, +1 because of the [) of range, total 0
                    try:
                        assignments[i] = str(r[2])
                    except IndexError:
                        self.log.error("residue %s not in sequence!" % i)
                        return None
            assignments = [ 'S' if a is None else a for a in assignments ]

        else:
            self.log.error("No disorder data or prediction available!")
            return None
        
        return assignments

    def _get_mobidb(self, sequence, use_alias=None):
        if use_alias is not None:
            gene_id = sequence.aliases[use_alias]
            self.log.info("using alias %s as gene name" % sequence.aliases[use_alias])
        else:
            gene_id = sequence.gene_id

        try:
            url = ('/').join([self._requests_url, gene_id, 'consensus'])
            self.log.debug("fetching %s" % url)
            req = rq.get(url)
        except:
            self.log.error("Couldn't get data for %s" % gene_id)
            return None

        if req.status_code == 200:
            try:
                data = req.json()
            except:
                self.log.error("Couldn't parse data for %s" % gene_id)
                return None
        else:
            self.log.warning("Data requests didn't complete correctly for %s" % gene_id)
            return None

        return data

class ManualAnnotation(StaticSource):

    _expected_cols = ['name', 'site', 'type', 'function', 'reference']
    _ptm_keywords = ['ptm_cleavage', 'ptm_phosphorylation', 'ptm_ubiquitination', 'ptm_acetylation', 'ptm_sumoylation', 'ptm_nitrosylation', 'ptm_methylation']
    _supported_position_properties = _ptm_keywords
    _supported_sequence_properties = ['linear_motif']

    @logger_init
    def __init__(self, datafile, **parsing_options):

        description="Annotations from %s" % datafile
        super(ManualAnnotation, self).__init__(name='Manual annotations', version='', description=description)

        self._datafile = datafile
        self._df = None
        
        try:
            self._parse_datafile(**parsing_options)
        except:
            self.log.error("An error occurred while parsing the input file %s; manual annotation will be skipped" % self._datafile)
            self._df = None
            return

    def _parse_datafile(self, **parsing_options):

        self.log.info("Parsing annotations from %s" % self._datafile)

        try:
            df = pd.read_csv(self._datafile, **parsing_options)
        except IOError:
            self.log.error("Parsing of file %s failed. ")
            raise IOError

        try:
            self._df = df[ self._expected_cols ]
        except IndexError:
            log.error("required columns not found in csv file (these are: %s). Manual annotation will be skipped" % ", ".join(self._expected_cols))
            raise IndexError

        all_properties = set(self._supported_position_properties + self._supported_sequence_properties)

        diff = set(df['type']).difference(all_properties)
        if len(diff) > 0:
            self.log.warning("the following annotation types were not recognized: %s" % ", ".join(diff))

    def add_mutations(self, sequence):

        if self._df is None:
            return

        tmp_df = self._df[ self._df['type'] == 'mutation' ]

        unique_mutations = list(set(tmp_df['site']))

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
                self.log.warning("for mutation %s, residue %d is %s in wild-type sequence; it will be skipped" %(m, num, position.wt_residue_type))
                continue

            mutation_obj = Mutation(sequence.positions[site_seq_idx],
                                    mut,
                                    [self])

            position.add_mutation(mutation_obj)

    def add_position_properties(self, sequence, prop=_supported_position_properties):

        if self._df is None:
            return

        tmp_df = self._df[ self._df.apply(lambda x: x['type'] in self._supported_position_properties, axis=1) ]

        for idx, row in tmp_df.iterrows():

            this_position = sequence.positions[sequence.seq2index(int(row['site']))]
            
            property_obj = position_properties_classes[row['type']](sources=[self],
                                                                    position=this_position)

            property_obj.metadata['function'] = row['function']
            property_obj.metadata['reference']      = row['reference']

            this_position.add_property(property_obj)

            self.log.info("added %s from row %d" % (row['type'], idx+1))


    def add_sequence_properties(self, sequence, prop=_supported_sequence_properties):

        if self._df is None:
            return

        tmp_df = self._df[ self._df.apply(lambda x: x['type'] in self._supported_sequence_properties, axis=1) ]

        for idx,row in tmp_df.iterrows():
            if '-' in row['site']:
                tmp = list(map(int, row['site'].split("-")))
                try:
                    assert(len(tmp) == 2)
                    assert(int(tmp[0]) <= int(tmp[1]))
                except AssertionError:
                    log.error("line %d range in site column must be specified as 'from-to' (ex 10-15)" % (idx+1))
                positions = tuple(range(tmp[0], tmp[1]+1))
            else:
                positions = ( int(row['site']) )

            try: 
                positions = [ sequence.positions[sequence.seq2index(p)] for p in positions ]
            except:
                self.log.warning("position property %s is outside the protein sequence; it will be skipped" % row['type'])
                continue

            property_obj = sequence_properties_classes[row['type']](sources=[self],
                                                                    positions=positions,
                                                                    lmtype=row['name'])


            property_obj.metadata['function'] = row['function']
            property_obj.metadata['reference']      = row['reference']

            sequence.add_property(property_obj)
