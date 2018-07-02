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
    def __init__(self, *args, **kwargs):
        description = "Uniprot knowledge-base"
        super(UniProt, self).__init__(name='UniProt', version='1.0', description=description)
        self._uniprot_service = bsUniProt()

    def get_sequence(self, gene_id):
        upid = self._uniprot_service.search("gene:%s AND organism:human" % gene_id, columns="entry name").split()[2]
        sequence = self._uniprot_service.get_fasta_sequence(str(upid))
        return Sequence(gene_id, sequence, self)

class cBioPortal(DynamicSource, object):
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
        print "Unique mutations found in %s: %s" % (self.name, ", ".join(sorted(unique_mutations, key=lambda x: int(x[1:-1]))))
        for m_idx,m in enumerate(unique_mutations):
            wt = m[0]
            mut = m[-1]
            num = int(m[1:-1])
            if wt == mut:
                continue
            mutation_indices = [i for i,x in enumerate(mutations) if x == m]
            try:
                site_seq_idx = sequence.seq2index(num)
            except:
                print "WARNING: residue %d out of range %s" % (num, out_metadata)
                for mi in mutation_indices:
                    print "£", out_metadata['cancer_study'][mi]
                continue
            position = sequence.positions[site_seq_idx]
            if position.wt_residue_type != wt:
                print "WARNING: Skipping %s %s (it's %s in sequence!)" % (m, self.name, position.wt_residue_type)
                for mi in mutation_indices:
                    print "£", out_metadata['cancer_study'][mi]
                continue

            mutation_obj = Mutation(sequence.positions[site_seq_idx],
                                    mut,
                                    [self])
            for md in metadata:
                mutation_obj.metadata[md] = []
                for mi in mutation_indices:
                    tmp_md = [self] + out_metadata[md][mi]
                    this_md = metadata_classes[md](*tmp_md)
                    mutation_obj.metadata[md].append(this_md)

            sequence.positions[site_seq_idx].add_mutation(mutation_obj)

    def _get_cancer_types(self):
        response = rq.post(self._requests_url, {'cmd':'getTypesOfCancer'}).text
        tmp = response.split("\n")
        data = []
        for line in tmp[:-1]:
            data.append(tuple(map(str, line.strip().split("\t"))))
        self._cancer_types = dict(data)

    def _get_cancer_studies(self):
        if self._cache_cancer_studies is None:
            self._cache_cancer_studies = []
            response = rq.post(self._requests_url, {'cmd':'getCancerStudies'}).text
            tmp = response.split("\n")
            for line in tmp[1:]:
                if line:
                    tmp2 = line.split("\t")
                    self._cache_cancer_studies.append(tuple(tmp2))

    def _get_genetic_profiles(self, cancer_studies=None):
        if self._cache_cancer_studies is None:
            if cancer_studies is not None:
                self._cache_cancer_studies = cancer_studies
            else:
                self._get_cancer_studies()
        
        self._cache_genetic_profiles = {}
        for cancer_study_id in zip(*self._cache_cancer_studies)[0]:
            self._cache_genetic_profiles[cancer_study_id] = []
            response = rq.post(self._requests_url, {'cmd':'getGeneticProfiles',
                                                    'cancer_study_id':cancer_study_id}).text
            tmp = response.split("\n")
            for line in tmp[1:]:
                if line:
                    tmp2 = line.split("\t")
                    if tmp2[4] == 'MUTATION' or tmp2[4] == 'MUTATION_EXTENDED':
                        self._cache_genetic_profiles[cancer_study_id].append((tmp2[0],tmp2[4]))

    def _get_case_sets(self, cancer_studies=None):
        if self._cache_cancer_studies is not None:
            if cancer_studies is not None:
                self._cache_cancer_studies = cancer_studies
            else:
                self._get_cancer_studies()

        self._cache_case_sets = {}
        for cancer_study_id in zip(*self._cache_cancer_studies)[0]:
            self._cache_case_sets[cancer_study_id] = []
            response = rq.post(self._requests_url, {'cmd':'getCaseLists',
                                                    'cancer_study_id':cancer_study_id}).text
            tmp = response.split("\n")
            for line in tmp[1:]:
                if line:
                    tmp2 = line.split("\t")
                    self._cache_case_sets[cancer_study_id].append(tmp2[0])

    def _get_profile_data(self, gene_id, cancer_studies=None, metadata=[], mapisoform=None, mainisoform=1):

        out_metadata = dict(zip(metadata, [list() for i in range(len(metadata))]))

        if self._cache_profile_data is None:
            self._get_genetic_profiles(cancer_studies=cancer_studies)
        if self._cache_case_sets is None:
            self._get_case_sets(cancer_studies=cancer_studies)

        if len(self._cache_cancer_studies[0]) == 3:
            for i,c in enumerate(self._cache_cancer_studies):
                self._cache_cancer_studies[i] = c + tuple([mainisoform])    

        mut_regexp = '[A-Z][0-9]+[A-Z]$'
        mut_prog = re.compile(mut_regexp)
        split_regexp = '\s+|,'
        split_prog = re.compile(split_regexp)
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

        for cancer_study_id, short_desc, long_desc, cancer_study_isoform in self._cache_cancer_studies:
            print "doing ", cancer_study_id, time.strftime('%x %X')
            if do_cancer_type:
                cancer_study_suffix = cancer_study_id.split("_")[0]
                try:
                    cancer_type = self._cancer_types[cancer_study_suffix]
                except KeyError:
                    cancer_type = cancer_study_suffix
            for case_set in self._cache_case_sets[cancer_study_id]:
                print "    doing", case_set, time.strftime('%x %X')
                cancer_study_mutation_ids = [x[0] for x in self._cache_genetic_profiles[cancer_study_id] if x[1] == 'MUTATION']
                if len(cancer_study_mutation_ids) > 0:
                    response = rq.post(self._requests_url, {'cmd':'getProfileData',
                                                            'case_set_id':case_set,
                                                            'genetic_profile_id':",".join(cancer_study_mutation_ids),
                                                            'gene_list':gene_id}).text

                    for line in response.strip().split("\n"):
                        if not line.startswith("#") and not line.startswith("GENE_ID"):
                            tmp = split_prog.split(line)[2:]
                            tmp2 = filter(lambda x: x != 'NaN', tmp)
                            tmp3 = filter(lambda x: mut_prog.match(x), tmp2)
                            mutations.extend(tmp3)
                            if do_cancer_type:
                                out_metadata['cancer_type'].extend([[cancer_type]]*len(tmp3))
                            if do_cancer_study:
                                out_metadata['cancer_study'].extend([[cancer_study_id]]*len(tmp3))
                            if do_genomic_coordinates:
                                out_metadata['genomic_coordinates'].extend([None]*len(tmp3))


                cancer_study_mutation_ids = [x[0] for x in self._cache_genetic_profiles[cancer_study_id] if x[1] == 'MUTATION_EXTENDED']
                if len(cancer_study_mutation_ids) > 0:
                    response = rq.post(self._requests_url, {'cmd':'getMutationData',
                                                            'case_set_id':case_set,
                                                            'genetic_profile_id':",".join(cancer_study_mutation_ids),
                                                            'gene_list':gene_id}).text

                    for line in response.strip().split("\n"):
                        if not line.startswith("#") and not line.startswith("entrez_gene_id") and line:
                            tmp=line.strip().split("\t")

                            try:
                                if tmp[5] != "Missense_Mutation":
                                    continue
                            except:
                                pass
                            if not mut_prog.match(tmp[7]):
                                continue

                            if cancer_study_isoform != mainisoform:
                                resn = int(tmp[7][1:-1])
                                tmp[7] = "%s%d%s" % (tmp[7][0], 
                                                    mapisoform[cancer_study_isoform][resn],
                                                    tmp[7][-1])
                                print "        converting isoform, %d to %d" % (resn, mapisoform[cancer_study_isoform][resn])
                            mutations.append(tmp[7])
                            gd = ['hg19', tmp[12], tmp[13], tmp[14], tmp[15]]

                            if do_cancer_type:
                                out_metadata['cancer_type'].append([cancer_type])
                            if do_cancer_study:
                                out_metadata['cancer_study'].append([cancer_study_id])
                            if do_genomic_coordinates:
                                out_metadata['genomic_coordinates'].append(gd)

        return mutations, out_metadata

class COSMIC(DynamicSource, object):
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

        for f in self._database_files:
            with open(f) as fh:
                fh.next()
                for line in fh:
                    tmp = line.strip().split("\t")
                    if tmp[0] != gene_id or tmp[11] not in cancer_types:
                        continue

                    if filter_snps and tmp[25] == 'y':
                        continue

                    mut_str = tmp[18]
                    if self._mut_prog.match(mut_str):
                        mutations.append(mut_str)
                        if do_cancer_type:
                            out_metadata['cancer_type'].append([tmp[11]])
                        if do_genomic_coordinates:
                            gd = []
                            if tmp[22] == '38':
                                gd.append('hg38') # genome version
                            elif tmp[22] == '37':
                                gd.append('hg19')
                            else:
                                out_metadata['genomic_coordinates'].append(['', '', '', '', ''])
                                continue
                            tmp2 = tmp[23].split(":")
                            gd.append(tmp2[0]) # chr
                            gd.extend(tmp2[1].split("-")) # [start, end]
                            gd.append(tmp[17][-3]) # ref
                            out_metadata['genomic_coordinates'].append(gd)

        return mutations, out_metadata

    def add_mutations(self, sequence, cancer_types=None, use_alias=None, filter_snps=True, metadata=[]):

        if cancer_types is None:
            cancer_types = self._get_all_cancer_types()
        if use_alias is not None:
            gene_id = sequence.aliases[use_alias]
        else:
            gene_id = sequence.gene_id
        raw_mutations, out_metadata = self._parse_db_file(gene_id, cancer_types=cancer_types, metadata=metadata)
        mutations = map(lambda x: x[2:], raw_mutations)
        unique_mutations = list(set(mutations))

        print "Unique mutations found in %s: %s" % (self.name, ", ".join(sorted(unique_mutations, key=lambda x: int(x[1:-1]))))

        for m in unique_mutations:
            wt = m[0]
            mut = m[-1]
            num = int(m[1:-1])
            if wt == mut:
                continue
            site_seq_idx = sequence.seq2index(num)
            position = sequence.positions[site_seq_idx]
            if position.wt_residue_type != wt:
                print "WARNING: Skipping %s %s (it's %s in sequence!)" % (m, self.name, position.wt_residue_type)

            mutation_indices = [i for i, x in enumerate(mutations) if x == m]

            mutation_obj = Mutation(sequence.positions[site_seq_idx],
                                    mut,
                                    [self])

            print mutation_indices
            for md in metadata:
                mutation_obj.metadata[md] = []
                for mi in mutation_indices:
                    tmp_md = [self] + out_metadata[md][mi]
                    this_md = metadata_classes[md](*tmp_md)
                    mutation_obj.metadata[md].append(this_md)
            position.add_mutation(mutation_obj)

class PhosphoSite(DynamicSource, object):
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
                site_seq_idx = sequence.seq2index(site)
                position = sequence.positions[site_seq_idx]
                if position.wt_residue_type != wt:
                    print "Warning: %s is %s in sequence!)" % (m, position.wt_residue_type)

                already_annotated = False
                for prop in position.properties:
                    if isinstance(prop, position_properties_classes[ptm]):
                        prop.sources.append(self)
                        already_annotated = True
                
                if not already_annotated:
                    property_obj = position_properties_classes[ptm](  sources=[self],
                                                        position=sequence.positions[site_seq_idx]
                                                        )
                    position.add_property(property_obj)

class MyVariant(DynamicSource, object):
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
            return False

        mutation.metadata['revel_score'] = list()

        gcs = list(set(mutation.metadata['genomic_coordinates']))

        for gc in gcs:
            print "adding revel to %s" % mutation
            if gc.coord_start != gc.coord_end:
                continue
            if gc.genome_version == 'hg38':
                converted_coords = self._lo.convert_coordinate('chr%s' % gc.chr, int(gc.coord_start))
                assert len(converted_coords) == 1
                converted_coords = (converted_coords[0][0], converted_coords[0][1])
            elif gc.genome_version == 'hg19':
                converted_coords = ('chr%s' % gc.chr, int(gc.coord_start))
            else:
                return False

            query_str = 'chr%s:%d' % converted_coords
            query = self._mv.query(query_str)
            hits = query['hits']
            if len(hits) < 1:
                continue
            revel_score = None
            for hit in hits:
                try:
                    if not mutation.sequence_position.sequence_position == hit['dbnsfp']['aa']['pos'] and \
                        not mutation.sequence_position.wt_residue_type   == hit['dbnsfp']['aa']['ref']:
                        raise TypeError
                except KeyError:
                    continue
                except TypeError:
                    try:
                        if not mutation.sequence_position.sequence_position in map(int, hit['dbnsfp']['aa']['pos']) and \
                            not mutation.sequence_position.wt_residue_type   == hit['dbnsfp']['aa']['ref']:
                            raise TypeError
                    except KeyError:
                        continue
                    except TypeError:
                        continue
                
                revel_score = None
                if mutation.mutated_residue_type == hit['dbnsfp']['aa']['alt']:
                    try:
                        revel_score = hit['dbnsfp']['revel']['score']
                    except KeyError:
                        continue

                if revel_score:
                    break

            if not revel_score:
                return False

            revel_md = DbnsfpRevel(source=self, score=revel_score)
            mutation.metadata['revel_score'] = [revel_md]
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
    def __init__(self):
        description = "ELM Prediction"
        super(ELMPredictions, self).__init__(name='ELM', version='1.0', description=description)


        self._classes_url  = "http://elm.eu.org/elms/"
        self._requests_url = "http://elm.eu.org/start_search/"
        self._get_elm_classes()

    def _get_elm_classes(self):
        self._elm_classes = {}
        unquote = lambda x: str.strip(x, '"')

        response = rq.get(self._classes_url+'/elms_index.tsv').text
        tmp = response.split("\n")
        for line in tmp:
            if line.startswith('"ELME'):
                tmp2 = line.strip().split("\t")
                tmp2 = map(str, tmp2)
                tmp2 = map(unquote, tmp2)
                self._elm_classes[tmp2[1]] = tmp2[2:]

    def _get_prediction(self, gene_name):
        req_url = os.path.join(self._requests_url, gene_name) + ".tsv"

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

        if use_alias is None:
            data = self._get_prediction(sequence.gene_id)
        else:
            data = self._get_prediction(sequence.aliases[use_alias])

        for d in data:
            if d[5]:
                continue

            if re.match(exclude_elm_classes, d[0]):
                continue

            this_positions = []
            for p in range(d[1],d[2]+1):
                this_positions.append(sequence.positions[sequence.seq2index(p)])

            if True:
                property_obj = sequence_properties_classes['linear_motif']  (sources=[self],
                                                             positions=this_positions,
                                                             lmtype=self._elm_classes[d[0]][0])

            else:
                continue
            property_obj.metadata['function'] = [self._elm_classes[d[0]][0]]
            property_obj.metadata['ref']      = self.description
            sequence.add_property(property_obj)

class ManualAnnotation(StaticSource):
    def __init__(self, datafile, csvoptions=None):
        description="Annotations from %s" % datafile
        super(ManualAnnotation, self).__init__(name='Manual annotations', version='', description=description)

        self._datafile = datafile
        fh = open(self._datafile, 'rb')
        if csvoptions is None:
            csvoptions = {}
        self._csv = csv.reader(fh, **csvoptions)
        self._table = None

    def _parse_datafile(self):
        out = [] # type seq type function ref
        ptm_keywords = ['cleavage', 'phosphorylation', 'ubiquitination', 'acetylation', 'sumoylation', 's-nitrosylation', 'methylation']
        for idx,row in enumerate(self._csv):
            if len(row) < 5:
                continue
            if   row[2] in ptm_keywords:
                this_row = [row[2], tuple([int(row[1])]), row[0], row[3], row[4]]
            elif row[2] == 'linear_motif':
                this_row = [row[2]]
                if '-' in row[1]:
                    tmp = map(int, row[1].split("-"))
                    this_row.append(tuple(range(tmp[0], tmp[1]+1)))
                else:
                    this_row.append( (int(row[1]) ))
                this_row.extend([row[0], row[3], row[4]])
            elif row[2] == 'mutation':
                this_row = [ row[2], tuple([int(row[1][1:-1])]), row[1], row[0], row[3], row[4] ]
            else:
                continue

            out.append(tuple(this_row))
        self._table = out


    def add_mutations(self, sequence):
        if self._table is None:
            self._parse_datafile()

        mutations = [d[2] for d in self._table if d[0] == 'mutation']
        unique_mutations = list(set(mutations))

        print "Unique mutations found in %s: %s" % (self.name, ", ".join(sorted(unique_mutations, key=lambda x: int(x[1:-1]))))

        for m in unique_mutations:
            wt = m[0]
            mut = m[-1]
            num = int(m[1:-1])
            if wt == mut:
                continue
            site_seq_idx = sequence.seq2index(num)
            position = sequence.positions[site_seq_idx]
            if position.wt_residue_type != wt:
                print "WARNING: Skipping %s %s (it's %s in sequence!)" % (m, self.name, position.wt_residue_type)
                continue

            mutation_indices = [i for i, x in enumerate(mutations) if x == m]

            mutation_obj = Mutation(sequence.positions[site_seq_idx],
                                    mut,
                                    [self])

            position.add_mutation(mutation_obj)

    def add_position_properties(self, sequence):
        if self._table is None:
            self._parse_datafile()

        for d in self._table:
            this_position = sequence.positions[sequence.seq2index(d[1][0])]
            
            try:
                property_obj = position_properties_classes[d[0]](sources=[self],
                                                                 position=this_position)
                assert len(d[1]) == 1
            except KeyError:
                continue

            property_obj.metadata['function'] = [d[3]]
            property_obj.metadata['ref']      = [d[4]]
            this_position.add_property(property_obj)

    def add_sequence_properties(self, sequence):
        if self._table is None:
            self._parse_datafile()

        for d in self._table:
            this_positions = []
            for p in d[1]:
                this_positions.append(sequence.positions[sequence.seq2index(p)])

            try:
                property_obj = sequence_properties_classes[d[0]](sources=[self],
                                                             positions=this_positions,
                                                             lmtype=d[2])
            except KeyError:
                continue
            property_obj.metadata['function'] = [d[3]]
            property_obj.metadata['ref']      = [d[4]]
            sequence.add_property(property_obj)



