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
from Bio.PDB.Polypeptide import three_to_one
from Bio import SeqIO
import numpy as np
import pandas as pd
from .core import Sequence, Mutation
from .properties import *
from .metadata import *
from .log import *
from io import StringIO
import json
import re
import sys
import os
import myvariant
import pyliftover
import itertools
from parse import parse
from multiprocessing.dummy import Pool as threadPool
from future.utils import iteritems
from bravado.client import SwaggerClient
from bravado.requests_client import RequestsClient
from requests.adapters import HTTPAdapter


import sys
if sys.version_info[0] >= 3:
    unicode = str

class Source(object):
    """Base class for implementing data sources.

    Attributes
    ----------

    name : :obj:`str`
        name of the data source
    version : :obj:`str`
        version of the data source
    description : :obj:`str`
        short description of the data source
    """
    def __init__(self, name, version, description):
        """Class constructor

        Parameters
        ----------
        name : :obj:`str`
            name of the data source
        version : :obj:`str`
            version of the data source
        description : :obj:`str`
            short description of the data source
        """

        self.name = name
        self.version = version
        self.description = description
    def get_sequence(self, gene_id):
        return None
    def get_mutations(self, gene_id):
        return None
    def get_position_properties(self, position):
        return None
    def get_mutation_properties(self, mutation):
        return None

class DynamicSource(Source, object):
    """Base class for implementing dynamic data sources. Dynamic data sources
    are remote data sources that can be queried through the internet."""

    def __init__(self, *args, **kwargs):
        super(DynamicSource, self).__init__(*args, **kwargs)
class StaticSource(Source, object):
    """Base class for implementing static data sources. Static data sources
    are local to the system in use and need to be provided manually. They
    don't change without manual intervention"""

    def __init__(self, *args, **kwargs):
        super(StaticSource, self).__init__(*args, **kwargs)

class UniProt(DynamicSource, object):
    """Class for the UniProt data source. It is used to download the protein
    sequence of the main UniProt isoform for a certain gene and build the 
    Sequence object, which is the main entry point for annotations in cancermuts"""

    @logger_init
    def __init__(self, *args, **kwargs):
        description = "Uniprot knowledge-base"
        super(UniProt, self).__init__(name='UniProt', version='1.0', description=description)
        self._uniprot_service = bsUniProt()

    def get_sequence(self, gene_id, upid=None):
        if upid is None:
            self.log.info("retrieving UniProt ID for human gene %s" % gene_id)
            try:
                uniprot_table = pd.read_csv(StringIO(self._uniprot_service.search(f"(gene:{gene_id}) AND (organism_id:9606)")), sep='\t')
                upids = uniprot_table['Entry Name'].to_list()
            except:
                self.log.error("Failed to retrieve list of Uniprot IDs")
                return None

            upid = upids[0]
            if len(upids) > 1:
                self.log.warning("the following UniProt entries were found for gene %s: %s; will use %s" %(gene_id, ', '.join(upids), upid))
            else:
                self.log.info("will use Uniprot ID %s" % upid)
        else:
            self.log.info("The user-provided Uniprot ID (%s) will be used" % upid)
        aliases = {'uniprot' : upid,
                   'entrez' :       self._get_aliases(upid, ['GeneID'])['GeneID'],
                   'uniprot_acc' : self._get_aliases(upid, ['UniProtKB_primaryAccession'])['UniProtKB_primaryAccession']}

        self.log.info("retrieving sequence for UniProt sequence for Uniprot ID %s, gene %s" % (upid, gene_id))

        try:
            sequence = self._get_fasta_sequence(str(upid))
        except:
            self.log.error("failed retrieving sequence for Uniprot ID %s" % upid)
            return None

        return Sequence(gene_id, sequence, self, aliases=aliases)

    def _get_fasta_sequence(self, upid):
        """This function downloads the sequence of the specified UniProt entry
        by means of the UniProt API. The downloaded sequence corresponds to the
        sequence of the main isoform.

        Parameters
        ----------
        upid : :obj:`str`
            a valid UniProt identifier

        Returns
        -------
        sequence : :obj:`str`
            a protein sequence as a string
        """

        try:
            fasta_text = self._uniprot_service.get_fasta(upid)
        except:
            self.log.error(f"could not retrieve FASTA sequence for {upid}")
            return None

        try:
            sequences = list(SeqIO.parse(StringIO(fasta_text), 'fasta'))
        except:
            self.log.error(f"could not parse obtained FASTA file for {upid}")
            return None

        if len(sequences) > 1:
            self.log.warning(f"Multiple FASTA sequences for {upid} found; the first one will be used")

        return str(sequences[0].seq)

    def _get_aliases(self, gene_id, to, fr='UniProtKB_AC-ID'):

        out = {}

        for t in to:

            if re.match('UniProtKB_\S+', t):
                t_keyword = "_".join(t.split("_")[1:])
                t_query = 'UniProtKB'
            else:
                t_keyword = None
                t_query = t

            responses = self._uniprot_service.mapping( fr = fr,
                                                       to = t_query,
                                                       query = gene_id )['results']

            if len(responses) == 0:
                self.log.warning(f"No {t} found for {gene_id}")
                return None
            if len(responses) > 1:
                raise TypeError
            if len(responses) > 1:
                self.log.warning(f"More than one {t} found: {' '.join(responses)}")
                self.log.warning("The first one will be used")

            results = responses[0]

            if t_keyword is not None:
                out[t] = results['to'][t_keyword]
            else:
                out[t] = results['to']
        return out


class cBioPortal(DynamicSource, object):

    default_strand = '+'
    
    @logger_init
    def __init__(self, cancer_studies=None):
        description = "cBioPortal"

        super(cBioPortal, self).__init__(name='cBioPortal', version='1.0', description=description)

        self._api_endpoint = 'https://www.cbioportal.org/api/v2/api-docs'

        try:    # no validation as recommended on the cBioPortal website
            self._client = SwaggerClient.from_url(self._api_endpoint,
                                                config={"validate_requests":False,
                                                        "validate_responses":False,
                                                        "validate_swagger_spec": False})
        except Exception as e:
            self.log.error("Could not initialize the Swagger client for cBioPortal")
            return None

        self._cache_cancer_studies = pd.DataFrame()
        self._cache_genetic_profiles = pd.DataFrame()
        self._cache_case_sets = pd.DataFrame()
        self._cache_profile_data = pd.DataFrame()
        self._cancer_types = pd.DataFrame()
        self._get_cancer_types()
        self._get_cancer_studies(cancer_studies)
        self._get_molecular_profiles()

        self._mut_regexp = '[A-Z][0-9]+[A-Z]$'
        self._mut_prog = re.compile(self._mut_regexp)

    @property
    def cancer_types(self):
        return self._cancer_types

    @property
    def cancer_studies(self):
        return self._cache_cancer_studies

    def add_mutations(self, sequence, metadata=[]):
        mutations, out_metadata = self._get_available_mutations(sequence.aliases['entrez'], metadata=metadata)
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
                self.log.warning(f"mutation {m} is outside the protein sequence; it will be skipped")
                continue

            position = sequence.positions[site_seq_idx]

            if position.wt_residue_type != wt:
                self.log.warning(f"for mutation {m}, residue {wt} is {position.wt_residue_type} in wild-type sequence; it will be skipped")
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
            result = self._client.Cancer_Types.getAllCancerTypesUsingGET().result()
        except:
            self.log.error("failed retrieving cancer types")
            self._cancer_types = None

        self._cancer_types = pd.DataFrame(dict(
            [ (attr, [ getattr(entry, attr) for entry in result ]) for attr in dir(result[0]) ]))
 
    def _get_cancer_studies(self, cancer_studies=None):
        self.log.info("retrieving cancer studies...")
        try:
            result = self._client.Studies.getAllStudiesUsingGET().result()
        except:
            self.log.error("failed retrieving cancer studies")
            self._cache_cancer_studies = None
            return

        df = pd.DataFrame(dict(
            [ (attr, [ getattr(entry, attr) for entry in result ]) for attr in dir(result[0]) ]))

        if cancer_studies is None:
            self.log.info("all cancer studies will be considered")
            self._cache_cancer_studies = df
        else:
            set_cs = set(cancer_studies)
            set_csid = set(df.studyId)

            if not set_cs.issubset(set_csid):
                notcommon = set_cs - set_csid
                self.log.warning("the following cancer studies are not in the cBioPortal cancer studies list and will be skipped: %s" % (", ".join(sorted(list(notcommon)))))

            cancer_studies = set_cs & set_csid
            self._cache_cancer_studies = df[ df.apply(lambda x: x.studyId in cancer_studies, axis=1) ]

    def _get_molecular_profiles(self):
        self.log.info("retrieving molecular profiles")

        if self._cache_cancer_studies.empty or self._cache_cancer_studies is None:
            self.log.warning("no cancer studies are present, therefore molecular profles won't be downloaded")
            return

        molecularProfileFilter = { "studyIds" : self._cache_cancer_studies['studyId'].values.tolist()}

        result = self._client.Molecular_Profiles.fetchMolecularProfilesUsingPOST(molecularProfileFilter = molecularProfileFilter, 
                                                                                 projection = 'DETAILED').result()

        result_df = pd.DataFrame(dict(
            [ (attr, [ getattr(entry, attr) for entry in result ]) for attr in dir(result[0]) ]))

        result_df = result_df[(result_df.molecularAlterationType == 'MUTATION') | (result_df.molecularAlterationType == 'MUTATION_EXTENDED')]

        self._cache_genetic_profiles = result_df

    def _get_mutation_data(self, gene_id):

        mutationMultipleStudyFilter = { "entrezGeneIds": [ gene_id ],
                                        "molecularProfileIds": self._cache_genetic_profiles.molecularProfileId.values.tolist()
                                      }

        result = self._client.Mutations.fetchMutationsInMultipleMolecularProfilesUsingPOST(mutationMultipleStudyFilter=mutationMultipleStudyFilter, 
                                                                                             projection='DETAILED').result()

        df = pd.DataFrame(dict(
            [ (attr, [ getattr(entry, attr) for entry in result ]) for attr in dir(result[0]) ]))

        df = df[df['mutationType'] == 'Missense_Mutation']
        df = df[df['proteinChange'].str.match(self._mut_regexp, case=True)]

        return df

    def _get_available_mutations(self, gene_id, metadata=[]):

        mutation_data = self._get_mutation_data(gene_id)

        out_metadata = dict(list(zip(metadata, [list() for i in range(len(metadata))])))

        mutations = []

        if 'genomic_mutations' in metadata and not 'genomic_coordinates' in metadata:
            metadata.append('genomic_coordinates')

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

        for cancer_study in self._cache_cancer_studies.iterrows():

            cancer_study = cancer_study[1]

            self.log.debug(f"Gathering mutations for {cancer_study.studyId}")

            cancer_study_id = cancer_study.studyId
            cancer_type     = cancer_study.cancerType
            cancer_type_id  = cancer_study.cancerTypeId

            if do_cancer_type:
                if cancer_type_id is None:
                    self.log.warning("cancer type ID not found - cancer type ID will be inferred from study ID")
                    cancer_type_id = cancer_study_id.split('_')[0]

                if cancer_type is None:
                    try:
                        cancer_type = self._cancer_types[ self._cancer_types['cancerTypeId'] == cancer_type_id ]
                        if len(cancer_type) != 1:
                            raise KeyError
                        cancer_type = cancer_type.values[0][2]
                    except KeyError:
                        cancer_type = cancer_type_id
                        self.log.warning("cancer type for %s not found - cancer type ID will be used instead" % cancer_type_id)

            molecular_profile_ids = self._cache_genetic_profiles[self._cache_genetic_profiles.studyId == cancer_study.studyId]
            self.log.debug(f"Molecular profile ids: {molecular_profile_ids.molecularProfileId.tolist()}")

            cancer_study_mutation_ids     = molecular_profile_ids[ molecular_profile_ids['molecularAlterationType'] == 'MUTATION' ]
            cancer_study_mutation_ext_ids = molecular_profile_ids[ molecular_profile_ids['molecularAlterationType'] == 'MUTATION_EXTENDED' ]

            self.log.debug(f"Molecular profile ids with MUTATION: {cancer_study_mutation_ids.molecularProfileId.tolist()}")
            self.log.debug(f"Molecular profile ids with MUTATION_EXTENDED: {cancer_study_mutation_ext_ids.molecularProfileId.tolist()}")

            if len(cancer_study_mutation_ids) > 0:

                mutations_df = mutation_data[ mutation_data.molecularProfileId.isin(cancer_study_mutation_ids.molecularProfileId) ]
                self.log.debug(f"Found {len(mutations_df)} mutations for MUTATION molecular profiles")

                this_mutations = mutations_df['proteinChange'].tolist()

                mutations.extend(this_mutations)
                if do_cancer_type:
                    out_metadata['cancer_type'].extend([[cancer_type]]*len(this_mutations))
                if do_cancer_study:
                    out_metadata['cancer_study'].extend([[cancer_study_id]]*len(this_mutations))
                if do_genomic_coordinates:
                    out_metadata['genomic_coordinates'].extend([[None]]*len(this_mutations))
                if do_genomic_mutations:
                    out_metadata['genomic_mutations'].extend([[None]]*len(this_mutations))

            if len(cancer_study_mutation_ext_ids) > 0:

                mutations_df = mutation_data[ mutation_data.molecularProfileId.isin(cancer_study_mutation_ext_ids.molecularProfileId) ]

                self.log.debug(f"Found {len(mutations_df)} mutations for MUTATION EXTENDED molecular profiles")

                for row in mutations_df.iterrows():
                    row = row[1]
                    mutations.append(row['proteinChange'])

                    if do_cancer_type:
                        out_metadata['cancer_type'].append([cancer_type])
                    if do_cancer_study:
                        out_metadata['cancer_study'].append([cancer_study_id])
                    if do_genomic_coordinates or do_genomic_mutations:
                        gd = list(row[['chr',
                                    'startPosition',
                                    'endPosition',
                                    'referenceAllele']].values)

                        gd = ['hg19', str(gd[0]), str(int(gd[1])), str(int(gd[2])), str(gd[3])]
                        out_metadata['genomic_coordinates'].append(gd)
                    if do_genomic_mutations:
                        if row['startPosition'] == row['endPosition']:
                            gm_fmt = f"{row['chr']}:g.{row['startPosition']}{row['referenceAllele']}>{row['variantAllele']}"
                            gm = ['hg19', gm_fmt]
                        elif (len(row['referenceAllele']) == len(row['variantAllele'])) and \
                                ((row['endPosition'] - row['startPosition'] + 1) == len(row['variantAllele'])):
                            gm_fmt = f"{row['chr']}:g.{row['startPosition']}_{row['endPosition']}delins{row['variantAllele']}"
                            gm = ['hg19', gm_fmt]
                            #self.log.info("Added delins! " + gm_fmt)
                        else:
                            #self.log.warning("mutation corresponds to neither a substitution or delins, genomic mutation won't be annotated")
                            gm = None

                        out_metadata['genomic_mutations'].append(gm)

        return mutations, out_metadata

class COSMIC(DynamicSource, object):
    @logger_init
    def __init__(self, database_files=None, database_encoding=None):
        description = "COSMIC Database"
        super(COSMIC, self).__init__(name='COSMIC', version='v87', description=description)

        self._mut_regexp = 'p\.[A-Z][0-9]+[A-Z]$'
        self._mut_prog = re.compile(self._mut_regexp)
        self._mut_snv_regexp = '^[0-9]+:g\.[0-9+][ACTG]>[ACTG]'
        self._mut_snv_prog = re.compile(self._mut_snv_regexp)
        self._site_kwd = ['Primary site', 'Site subtype 1', 'Site subtype 2', 'Site subtype 3']
        self._histology_kwd = ['Primary histology', 'Histology subtype 1', 'Histology subtype 2', 'Histology subtype 3']
        self._use_cols_snp = [  'Gene name',
                                 'Mutation AA',
                                 'GRCh',
                                 'Mutation genome position',
                                 'Mutation CDS',
                                 'Mutation strand',
                                 'SNP',
                                 'HGVSG' ] + self._site_kwd + self._histology_kwd

        self._use_cols = [  'Gene name',
                            'Mutation AA',
                            'GRCh',
                            'Mutation genome position',
                            'Mutation CDS',
                            'Mutation strand',
                            'HGVSG' ] + self._site_kwd + self._histology_kwd

        dataframes = []

        if database_files is None:
            database_dir   = '/data/databases/cosmic'
            databases      = [  'CosmicMutantExport' ]
            self._database_files = [os.path.join(database_dir, i)+'.tsv' for i in databases]
        else:
            self._database_files = database_files

        if database_encoding is None or isinstance(database_encoding, str):
            encodings = [ database_encoding for i in self._database_files ]
        else:
            if len(database_encoding) != len(self._database_files):
                raise TypeError("encoding for COSMIC database files must be None, a single string, or "
                        "a list of strings, one per file")
            encodings = database_encoding

        for fi, f in enumerate(self._database_files):
            self.log.info("Parsing database file %s ..." % f)
            try:
                dataframes.append( pd.read_csv(f, sep='\t', dtype='str', na_values='NS', usecols=self._use_cols_snp, encoding=encodings[fi]) )
            except ValueError:
                try:
                    dataframes.append( pd.read_csv(f, sep='\t', dtype='str', na_values='NS', usecols=self._use_cols, encoding=encodings[fi]) )
                except ValueError:
                    self.log.error("Couldn't parse database file {}".format(fi))
                    self._df = None
                    return

        self._df = pd.concat(dataframes, ignore_index=True, sort=False)

    def _parse_db_file(self, gene_id, cancer_types=None,
                       cancer_histology_subtype_1=None,
                       cancer_histology_subtype_2=None,
                       cancer_histology_subtype_3=None,
                       cancer_sites=None,
                       cancer_site_subtype_1=None,
                       cancer_site_subtype_2=None,
                       cancer_site_subtype_3=None,
                       metadata=[], filter_snps=True):

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
            if 'SNP' in df.columns:
                df = df[ df['SNP'] != 'y' ]
            else:
                self.log.warning("Could not filter COSMIC for SNP as SNP column is not present")

        if cancer_types is not None:
            df = df[ df['Primary histology'  ].isin(cancer_types) ]
        if cancer_histology_subtype_1 is not None:
            df = df[ df['Histology subtype 1'].isin(cancer_histology_subtype_1) ]
        if cancer_histology_subtype_2 is not None:
            df = df[ df['Histology subtype 2'].isin(cancer_histology_subtype_2) ]
        if cancer_histology_subtype_3 is not None:
            df = df[ df['Histology subtype 3'].isin(cancer_histology_subtype_3) ]

        if cancer_sites is not None:
            df = df[ df['Primary site'  ].isin(cancer_sites) ]
        if cancer_site_subtype_1 is not None:
            df = df[ df['Site subtype 1'].isin(cancer_site_subtype_1) ]
        if cancer_site_subtype_2 is not None:
            df = df[ df['Site subtype 2'].isin(cancer_site_subtype_2) ]
        if cancer_site_subtype_3 is not None:
            df = df[ df['Site subtype 3'].isin(cancer_site_subtype_3) ]

        df = df[ df['Mutation AA'].notna() ]

        df = df[ df.apply(lambda x: bool(self._mut_prog.match(x['Mutation AA'])), axis=1) ]

        for r in df.iterrows():
            r = r[1]
            mutations.append(r['Mutation AA'])

            if do_cancer_type:
                out_metadata['cancer_type'].append([r['Primary histology']])

            if do_genomic_coordinates or do_genomic_mutations:
                gd = []
                grch = str(r['GRCh'])
                if grch == '38':
                    gd.append('hg38') # genome version
                elif grch == '37' or grch == '19':
                    gd.append('hg19')
                else:
                    gd = None
                    self.log.warning("Genome assembly not specified for mutation; genomic coordinates won't be annotated")

                if gd is not None:
                    tmp2 = r['Mutation genome position'].split(":")
                    gd.append(tmp2[0]) # chr
                    gd.extend(tmp2[1].split("-")) # [start, end]
                    gd.append(r['Mutation CDS'][-3]) # ref

            if do_genomic_coordinates:
                out_metadata['genomic_coordinates'].append(gd)

            if do_genomic_mutations:
                if gd is None:
                    self.log.warning("couldn't annotate genomic mutation")
                    gm = None
                else:
                    gm = [gd[0], r['HGVSG']]
                out_metadata['genomic_mutations'].append(gm)

            if do_site:
                out_metadata['cancer_site'].append([r.__getitem__(a) for a in self._site_kwd])

            if do_histology:
                out_metadata['cancer_histology'].append([r.__getitem__(a) for a in self._histology_kwd])

        return mutations, out_metadata

    def add_mutations(self, sequence, cancer_types=None,
                    cancer_histology_subtype_1=None,
                    cancer_histology_subtype_2=None,
                    cancer_histology_subtype_3=None,
                    cancer_sites=None,
                    cancer_site_subtype_1=None,
                    cancer_site_subtype_2=None,
                    cancer_site_subtype_3=None,
                    use_alias=None, filter_snps=True, metadata=[]):

        if cancer_types is None:
            self.log.info("no cancer type specified; will use all of them")
        if use_alias is not None:
            gene_id = sequence.aliases[use_alias]
            self.log.info("using alias %s as gene name" % sequence.aliases[use_alias])
        else:
            gene_id = sequence.gene_id

        raw_mutations, out_metadata = self._parse_db_file(gene_id, cancer_types=cancer_types,
                                                            cancer_histology_subtype_1=cancer_histology_subtype_1,
                                                            cancer_histology_subtype_2=cancer_histology_subtype_2,
                                                            cancer_histology_subtype_3=cancer_histology_subtype_3,
                                                            cancer_sites=cancer_sites,
                                                            cancer_site_subtype_1=cancer_site_subtype_1,
                                                            cancer_site_subtype_2=cancer_site_subtype_2,
                                                            cancer_site_subtype_3=cancer_site_subtype_3,
                                                            metadata=metadata)

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
                self.log.warning(f"mutation {m} is outside the protein sequence; it will be skipped")
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
    def __init__(self, database_dir, database_files=None):
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

        self._database_dir = database_dir

        if database_files is None:
            self._database_files = dict([(i, "%s/%s_site_dataset"%(self._database_dir, i.capitalize()[0] + i[1:])) for i in self._ptm_types])

        else:
            self._database_files = database_files

        self._dataframes = {}
        for k,f in iteritems(self._database_files):
            try:
                self._dataframes[k] = pd.read_csv(f, skiprows=3, sep='\t')
                self._dataframes[k] = self._dataframes[k][ self._dataframes[k]['ORGANISM'] == 'human']
            except:
                self.log.error("couldn't read database file %s" % f)
                raise IOError

    def _parse_db_file(self, gene_id):

        sites = dict(list(zip(self._ptm_types, [list() for i in range(len(self._ptm_types))])))

        for ptm_idx,ptm in enumerate(self._ptm_types):
            p_regexp = '[A-Z][0-9]+-%s' % self._ptm_suffixes[ptm_idx]
            p_prog = re.compile(p_regexp)
            p_sites = []
            df = self._dataframes[ptm]
            df = df[ df['GENE'].str.upper() == gene_id.upper() ]
            if df.empty:
                continue
            df = df[ df.apply(lambda x: bool(p_prog.match(x['MOD_RSD'])), axis=1) ]

            sites[ptm] = [ x[:self._ptm_suffix_offsets[ptm_idx]] for x in df['MOD_RSD'].values  ]

        return sites

    def add_position_properties(self, sequence, properties=None):

        if properties is None:
            properties = self._ptm_types

        sites = self._parse_db_file(sequence.gene_id)
        for ptm_idx,ptm in enumerate(properties):

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
                    self.log.warning("for PTM %s, residue %s is %s in wild-type sequence; it will be skipped" %(m, wt, position.wt_residue_type))
                    continue

                already_annotated = False
                for prop in position.properties:
                    if isinstance(prop, position_properties_classes[self._ptm_types_to_classes[ptm]]):
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
        self._revel_cache_gc = {}
        self._revel_cache_gm = {}


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
            self.log.debug(f"pulling revel score for {mutation}, {gc}, {gm}")
            if gm is not None:
                self.log.debug("genomic mutation will be used to retrieve revel score for mutation %s" % mutation)
                revel_score = self._get_revel_from_gm(mutation, gm)
            elif gc is not None:
                self.log.debug("genomic coordinates will be used to retrieve revel score for mutation %s" % mutation)
                revel_score = self._get_revel_from_gc(mutation, gc)
            else:
                self.log.warning("mutations %s has no genomic coordinates or genomic mutation available; it will be skipped" % mutation)

            self.log.debug(f"final revel score for {gc} {gm} {revel_score}")
            mutation.metadata['revel_score'].append(revel_score)

        self.log.debug(f"final revel metadata for {mutation} {mutation.metadata['revel_score']}")
    def _validate_revel_hit(self, mutation, hit):

        try:
            aa = hit['dbnsfp']['aa']
        except:
            self.log.error("no residue information was found; it will be skipped")
            return False

        if not isinstance(aa, list):
            aa = [aa]

        aa_short = []
        self.log.debug("{0} residue definitions will be tested".format(len(aa)))
        for idx_this_aa, this_aa in enumerate(aa):
            try:
                this_hit_pos = this_aa['pos']
                this_hit_ref = this_aa['ref']
                this_hit_alt = this_aa['alt']
            except:
                self.log.info("no residue information was found in variant definition {0}".format(idx_this_aa))
                continue

            if not mutation.sequence_position.wt_residue_type == this_hit_ref:
                self.log.info("reference residue in revel does not correspond in variant definition {0}".format(idx_this_aa))
                continue

            if not mutation.sequence_position.sequence_position == this_hit_pos:
                try:
                    if not mutation.sequence_position.sequence_position in [ int(a) for a in this_hit_pos ]:
                        raise TypeError
                except:
                    self.log.info("sequence position in revel does not correspond for in variant definition {0}".format(idx_this_aa))
                    continue

            if mutation.mutated_residue_type != this_hit_alt:
                self.log.info("protein mutation in revel does not correspond for this in variant definition {0}".format(idx_this_aa))
                continue

            aa_short.append(idx_this_aa)

            if len(aa_short) == 1:
                self.log.info("A single valid residue entry was found for revel score in mutation {0}".format(mutation))
                return True
            elif len(aa_short) > 1:
                self.log.warning("More than one valid residue entry found for revel score in mutation {0}".format(mutation))
                return True
            elif len(aa_short) == 0:
                self.log.warning("No valid residue entry found for revel score in mutation {0}".format(mutation))
            return False

    def _convert_hg38_to_hg19(self, gc):
        if gc.genome_build == 'hg38':
            self.log.info("%s has genomic data in hg38 assembly - will be converted to hg19" % gc)
            converted_coords = self._lo.convert_coordinate('chr%s' % gc.chr, int(gc.get_coord()))
            assert len(converted_coords) == 1
            converted_coords = (converted_coords[0][0], converted_coords[0][1])
        elif gc.genome_build == 'hg19':
            converted_coords = ('chr%s' % gc.chr, int(gc.get_coord()))
        else:
            self.log.error("genomic coordinates are not expressed either in hg38 or hg19 for %s; it will be skipped" % gc)
            converted_coords = None
        return converted_coords

    def _get_revel_from_gm(self, mutation, gm):

        if not gm.is_snv:
            return None

        converted_coords = self._convert_hg38_to_hg19(gm)
        if converted_coords is None:
            return None

        query_str = '%s:g.%d%s>%s' % (converted_coords[0], converted_coords[1], gm.ref, gm.alt)

        try:
            revel_score = self._revel_cache_gm[query_str]
            self.log.info(f"Revel score for {query_str} retrieved from cache")
            return DbnsfpRevel(source=self, score=revel_score)
        except KeyError:
            pass

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

        self._revel_cache_gm[query_str] = revel_score

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

    def _get_prediction(self, gene_name, elm_wait):
        
        if elm_wait > 0:
            self.log.info(f"waiting {elm_wait} seconds before querying ELM as requested")
        self.log.info("retrieving prediction for %s" % gene_name )
        time.sleep(elm_wait)
        
        try:
            req_url = os.path.join(self._requests_url, gene_name) + ".tsv"
            response = rq.get(req_url)

            if response.status_code == 429:
                self.log.error("connection to ELM server was refused as not enough time has passed since the last attempt; please wait at least 3 minutes before retrying")
                return None
            elif response.status_code != 200:
                self.log.error(f"ELM server responded with {response.status_code}; couldn't fetch prediction")
                return None
            else:
                response = response.text
        except:
            self.log.error("couldn't fetch ELM predictions")
            return None

        out = []
        
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

    def add_sequence_properties(self, sequence, exclude_elm_classes='{100}', use_alias='uniprot', elm_wait=180):
        self.log.info("adding ELM predictions to sequence ...")
        if use_alias is None:
            data = self._get_prediction(sequence.gene_id, elm_wait=elm_wait)
        else:
            self.log.info("will use alias %s as gene name" % sequence.aliases[use_alias])
            data = self._get_prediction(sequence.aliases[use_alias], elm_wait=elm_wait)

        for d in data:
            if d[5]:
                self.log.info("%s was filtered out by ELM" % d[0])
                continue

            if re.match(exclude_elm_classes, d[0]):
                self.log.info("%s was filtered out as requested" % d[0])
                continue

            this_positions = []
            for p in range(d[1],d[2]+1):
                this_positions.append(sequence.positions[sequence.seq2index(p)])

            property_obj = sequence_properties_classes['linear_motif']  (sources=[self],
                                                                         positions=this_positions,
                                                                         name=self._elm_classes[d[0]][0])

            property_obj.metadata['function'] = [self._elm_classes[d[0]][0]]
            property_obj.metadata['ref']      = self.description
            sequence.add_property(property_obj)

class gnomAD(DynamicSource, object):
    
    description = "gnomAD"

    _versions = {  '2.1' :             'gnomad_r2_1',
                   '3' :               'gnomad_r3',
                   '2.1_controls' :    'gnomad_2.1_controls',
                   '2.1_non-neuro' :   'gnomad_2.1_non_neuro',
                   '2.1_non-cancer' :  'gnomad_2.1_non_cancer',
                   '2.1_non-topmed' :  'gnomad_2.1_non_topmed',
                   'exac' :            'exac',
                }

    _assembly = {  '2.1' :             'GRCh37',
                   '3' :               'GRCh38',
                   '2.1_controls' :    'GRCh37',
                   '2.1_non-neuro' :   'GRCh37',
                   '2.1_non-cancer' :  'GRCh37',
                   '2.1_non-topmed' :  'GRCh37',
                   'exac' :            'GRCh37',
                }

    _version_str = {  '2.1' :              'gnomAD v2.1',
                       '3' :               'gnomAD v3',
                       '2.1_controls' :    'gnomAD v2.1 (controls)',
                       '2.1_non-neuro' :   'gnomAD v2.1 (non-neuro)',
                       '2.1_non-cancer' :  'gnomAD v2.1 (non-cancer)',
                       '2.1_non-topmed' :  'gnomAD v2.1 (non-topmed)',
                       'exac' :            'EXaC',
                    }


    @logger_init
    def __init__(self, version='2.1'):
        
        version = str(version)
        if version not in self._versions.keys():
            self.log.error("gnomAD version %s not supported by the current implementation" % version)
            raise TypeError

        super(gnomAD, self).__init__(name='gnomAD', version=version, description=self.description)

        self._gnomad_endpoint = 'https://gnomad.broadinstitute.org/api/'
        self._cache = {}
        self._supported_metadata = {'gnomad_exome_allele_frequency' : self._get_exome_allele_freq,
                                    'gnomad_genome_allele_frequency' : self._get_genome_allele_freq,
                                    'gnomad_popmax_exome_allele_frequency' : self._get_popmax_exome_allele_freq,
                                    'gnomad_popmax_genome_allele_frequency' : self._get_popmax_genome_allele_freq}

        gnomADExomeAlleleFrequency.set_version_in_desc(self._version_str[version])
        gnomADGenomeAlleleFrequency.set_version_in_desc(self._version_str[version])

    def add_metadata(self, sequence, md_type=['gnomad_exome_allele_frequency'], use_alias=None):
        
        if use_alias is not None:
            gene_id = sequence.aliases[use_alias]
            self.log.info("using alias %s as gene name" % sequence.aliases[use_alias])
        else:
            gene_id = sequence.gene_id

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
                self.log.warning("gnomAD doesn't support metadata type %s" % md_type)

        self.log.debug("collected metadata functions: %s" % ', '.join([i for i in metadata_functions.__repr__()]))

        for pos in sequence.positions:
            for mut in pos.mutations:
                for i,add_this_metadata in enumerate(metadata_functions):
                    mut.metadata[md_types[i]] = []
    
                if not 'genomic_mutations' in mut.metadata:
                    self.log.warning("no genomic mutation data available for gnomAD, mutation %s. It will be skipped" % mut)
                    continue
                elif mut.metadata['genomic_mutations'] is None or len(mut.metadata['genomic_mutations']) == 0:
                    self.log.warning("no genomic mutation data available for gnomAD, mutation %s. It will be skipped" % mut)
                    continue

                for i,add_this_metadata in enumerate(metadata_functions):
                    self.log.debug("adding metadata %s to %s" % (md_types[i], mut))
                    add_this_metadata(mut, gene_id)

    def _get_exome_allele_freq(self, mutation, gene_id):
        self.log.info("getting genome allele frequency")
        self._get_metadata(mutation, gene_id, 'gnomad_exome_allele_frequency')

    def _get_genome_allele_freq(self, mutation, gene_id):
        self.log.info("getting genome allele frequency")
        self._get_metadata(mutation, gene_id, 'gnomad_genome_allele_frequency')

    def _get_popmax_exome_allele_freq(self, mutation, gene_id):
        self.log.info("getting popmax genome allele frequency")
        self._get_metadata(mutation, gene_id, 'gnomad_popmax_exome_allele_frequency')

    def _get_popmax_genome_allele_freq(self, mutation, gene_id):
        self.log.info("getting popmax genome allele frequency")
        self._get_metadata(mutation, gene_id, 'gnomad_popmax_genome_allele_frequency')

    def _edit_variant_id(self, row):
        split_id = row['variant_id'].split('-')
        if len(split_id[2]) > 1 or len(split_id[3]) > 1:
            split_id[2] = '?'
        return "-".join(split_id)

    def _get_metadata(self, mutation, gene_id, md_type):

        exac_key = {    'gnomad_exome_allele_frequency'  : 'exome_af', 
                        'gnomad_genome_allele_frequency' : 'genome_af',
                        'gnomad_popmax_exome_allele_frequency' : 'popmax_exome_af',
                        'gnomad_popmax_genome_allele_frequency' : 'popmax_genome_af'       }

        ref_assembly = self._assembly[self.version]

        mutation.metadata[md_type] = list()

        if gene_id not in self._cache.keys():

            data = self._get_gnomad_data(gene_id, self._assembly[self.version], self._versions[self.version])

            if data is None:
                self._cache[gene_id] = None
                return None
            
            assemblies = set(data['reference_genome'])
            
            if len(assemblies) != 1:
                self.log.error("the downloaded data refers to more than one assembly!")
                return None

            if list(assemblies)[0] != ref_assembly:
                self.log.error("the downloaded data refers to an unexpected assembly!")
                return None

            self._cache[gene_id] = data
        
        else:
            data = self._cache[gene_id]

            if data is None:
                self.log.warning(f"cached data for gene {gene_id} was available, but didn't contain any information")
                return None

            self.log.info(f"data for gene {gene_id} already in cache")

        for variant in mutation.metadata['genomic_mutations']:
            if type(variant) is GenomicMutation and (variant.is_snv or variant.is_insdel):
                v_str = variant.as_assembly(ref_assembly).get_value_str(fmt='gnomad')
            else:
                v_str = None

            if variant.is_snv:
                this_df = data[ data['variant_id'] == v_str ]
            elif variant.is_insdel:
                this_df = data[ data['edited_variant_id'] == v_str ]
            else:
                self.log.info(f"variant not supported by gnomad {variant.description}")
                mutation.metadata[md_type].append(metadata_classes[md_type](self, None))
                continue

            if len(this_df) == 0:
                af = None
                self.log.info("no entry found for %s" % v_str)
            elif len(this_df) > 1:
                self.log.warning("more than one entry for %s! Skipping" % v_str)
                af = None
            elif len(this_df) == 1:
                if np.isnan(this_df[exac_key[md_type]].values[0]):
                    af = None
                    mutation.metadata[md_type] = np.nan
                    self.log.info("nan entry found for %s" % v_str)
                else:
                    af = this_df[exac_key[md_type]].values[0]
                    self.log.info("entry found for %s" % v_str)
            if af is not None:
                mutation.metadata[md_type].append(metadata_classes[md_type](self, af))

    def _get_gnomad_data(self, gene_id, reference_genome, dataset):
        headers = { "content-type": "application/json" }
        request="""query getGene($geneSymbol : String!,
                              $refBuild : ReferenceGenomeId!,
                              $dataset : DatasetId!) {
                gene(gene_symbol : $geneSymbol, reference_genome: $refBuild) {
                    gene_id
                    symbol
                    hgnc_id
                    variants(dataset : $dataset) {
                        variant_id
                        reference_genome
                        chrom
                        pos
                        ref
                        alt
                        rsids
                        exome {
                            ac
                            an
                            populations {
                                id
                                ac
                                an
                            }
                        }
                        genome {
                            ac
                            an
                            populations {
                                id
                                ac
                                an
                            }
                        }
                    }
                }
            }"""
        
        self.log.info("retrieving data for gene %s" % gene_id)

        try:
            response = rq.post(self._gnomad_endpoint,
                data=json.dumps({
                                "query": request,
                                "variables": { "geneSymbol": gene_id,
                                               "refBuild"  : reference_genome,
                                               "dataset"   : dataset }
                                }),
                headers={"Content-Type": "application/json"}) 
        except:
            self.log.error("Couldn't perform request for gnomAD")
            return None

        try:
            response_json = response.json()
        except:
            self.log.error("downloaded data couldn't be understood")
            return None

        if 'errors' in response_json.keys():
            self.log.error('The following errors were reported when querying gnomAD:')
            for e in response_json['errors']:
                self.log.error('\t%s' % e['message'])
            return None

        variants = pd.json_normalize(response_json['data']['gene']['variants'])
        variants = variants.rename(columns={'exome.ac':'exome_ac',
                                            'exome.an':'exome_an',
                                            'genome.ac':'genome_ac',
                                            'genome.an':'genome_an'
                                            })
        if variants.empty:
            self.log.warning('No variants were available for the specified gene')
            return None

        variants['edited_variant_id'] = variants.apply(self._edit_variant_id, axis=1)

        variants['total_ac'] = variants['exome_ac'] + variants['genome_ac']
        variants['total_an'] = variants['exome_an'] + variants['genome_an']

        variants['exome_af']  =  variants['exome_ac'] / variants['exome_an']
        variants['genome_af'] = variants['genome_ac'] / variants['genome_an']
        variants['total_af'] = variants['total_ac'] / variants['total_an']

        variants[['popmax_exome_ac',  'popmax_exome_an',  'popmax_exome_af',
              'popmax_genome_ac', 'popmax_genome_an', 'popmax_genome_af',
              'popmax_tot_ac',    'popmax_tot_an',    'popmax_tot_af']] = variants.apply(self.get_popmax_af, axis=1)

        return variants

    def get_popmax_af(self, variants, pops=['afr', 'eas', 'nfe', 'amr', 'sas']):
    # from https://gnomad.broadinstitute.org/help/popmax

    # This annotation contains allele frequency information (AC, AN, AF, homozygote
    # count) for the non-bottlenecked population with the highest frequency.

    # For gnomAD v2, this excludes
    #   Ashkenazi Jewish (asj),
    #   European Finnish (fin),
    #   and "Other" (oth) populations.
    # 
    # For gnomAD v3, this excludes
    #   Amish (ami),
    #   Ashkenazi Jewish (asj),
    #   European Finnish (fin),
    #   Middle Eastern (mid), and
    #   "Other" (oth) populations.

    # gnomAD version                  2    3
    ############################################################
    # afr: african/african american   +    +   
    # ami: Amish                      NA   -
    # amr: latino/admixed american    +    +
    # asj: ashkenazi jews             -    -
    # eas: east asian                 +    +
    # fin: european (finnish)         -    -
    # nfe: european (non-finnish)     +    +
    # mid: middle eastern             NA   -
    # oth: other                      -    -
    # sas: south asian                +    +
    
        popmax_allowed_pops = ['afr', 'amr', 'eas', 'nfe', 'sas']

        popmax_exome = pd.NA
        if 'exome.populations' in variants.keys().tolist():
            if type(variants['exome.populations']) == list:
                do_exome = True
                popmax_exome = pd.DataFrame(variants['exome.populations']).set_index('id').loc[popmax_allowed_pops,:]
                popmax_exome['af'] = popmax_exome['ac'] / popmax_exome['an']
                popmax_exome = popmax_exome.loc[popmax_exome['af'].idxmax()]

                popmax_exome_af = popmax_exome.af
                popmax_exome_ac = popmax_exome.ac
                popmax_exome_an = popmax_exome.an
            else:
                do_exome = False
        else:
            do_exome = False

        if do_exome == False:
            popmax_exome_ac = pd.NA
            popmax_exome_an = pd.NA
            popmax_exome_af = pd.NA

        
        popmax_genome = pd.NA
        if 'genome.populations' in variants.keys().tolist():
            if type(variants['genome.populations']) == list:
                do_genome = True
                popmax_genome = pd.DataFrame(variants['genome.populations']).set_index('id').loc[['afr', 'eas', 'nfe', 'amr', 'sas'],:]
                popmax_genome['af'] = popmax_genome['ac'] / popmax_genome['an']
                popmax_genome = popmax_genome.loc[popmax_genome['af'].idxmax()]

                popmax_genome_af = popmax_genome.af
                popmax_genome_ac = popmax_genome.ac
                popmax_genome_an = popmax_genome.an
    
            else:
                do_genome = False
        else:
            do_genome = False

        if do_genome == False:
                popmax_genome_ac = pd.NA
                popmax_genome_an = pd.NA
                popmax_genome_af = pd.NA


        if do_exome and do_genome:
            popmax_tot_ac = popmax_exome_ac + popmax_genome_ac
            popmax_tot_an = popmax_exome_an + popmax_genome_an
            popmax_tot_af = popmax_tot_ac / popmax_tot_an

        else:
            popmax_tot_ac = pd.NA
            popmax_tot_an = pd.NA
            popmax_tot_af = pd.NA

        return pd.Series([popmax_exome_ac,  popmax_exome_an,  popmax_exome_af,
                        popmax_genome_ac, popmax_genome_an, popmax_genome_af,
                        popmax_tot_ac,    popmax_tot_an,    popmax_tot_af]).replace(pd.NA,np.nan)



class MobiDB(DynamicSource):

    description = "MobiDB"

    @logger_init
    def __init__(self):

        super(MobiDB, self).__init__(name='MobiDB', version='0.1', description=self.description)

        self._requests_url = 'http://mobidb.bio.unipd.it/ws/'
        self._data_cache = {}
        self._supported_properties = { 'mobidb_disorder_propensity' : self._get_mobidb_disorder_predictions }
                                    
    def add_position_properties(self, sequence, prop=['mobidb_disorder_propensity'], use_alias='uniprot_acc'):
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

        if assignments is None:
            self.log.error("Couldn't retrieve prediction from MobiDB")
            return False

        if len(assignments) != len(sequence.positions): 
            self.log.error("Error in the assignment of predictions: length of MobiDB structure assignment and sequence don't match")
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
                self.log.warning("More than one prediction found for MobiDB; will use the first")

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

    def _get_mobidb(self, sequence, use_alias='uniprot_acc'):
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
            self.log.warning("Data requests for MobiDB didn't complete correctly for %s" % gene_id)
            return None

        return data

class ManualAnnotation(StaticSource):

    _expected_cols = ['name', 'site', 'type', 'function', 'reference']
    _ptm_keywords = ['ptm_cleavage', 'ptm_phosphorylation', 'ptm_ubiquitination', 'ptm_acetylation', 'ptm_sumoylation', 'ptm_nitrosylation', 'ptm_methylation']
    _supported_position_properties = _ptm_keywords
    _supported_sequence_properties = ['linear_motif', 'structure']
    _supported_mutation            = ['mutation']
    _mut_prot_parse = 'p.{wt:3l}{position:d}{mut:3l}'

    @logger_init
    def __init__(self, datafile, **parsing_options):

        description="Annotations from %s" % datafile
        super(ManualAnnotation, self).__init__(name=f"Manual annotations from {datafile}",
            version='',
            description=description)

        self._datafile = datafile
        self._df = None

        if "sep" not in parsing_options.keys():
            parsing_options['sep'] = ';'

        try:
            self._parse_datafile(**parsing_options)
        except pd.errors.ParserError:
            self.log.error("An error occurred while parsing the input file %s; manual annotation will be skipped" % self._datafile)
            self._df = None
            return

    def _parse_datafile(self, **parsing_options):

        self.log.info("Parsing annotations from %s" % self._datafile)

        try:
            df = pd.read_csv(self._datafile, **parsing_options, keep_default_na=False)
        except IOError:
            self.log.error("Parsing of file %s failed. ")
        try:
            self._df = df[ self._expected_cols ]
        except:
            self.log.error("required columns not found in csv file (these are: %s). Manual annotation will be skipped" % ", ".join(self._expected_cols))
            raise IndexError

        self.log.info("Parsed file:")
        self.log.info('\n{0}'.format(str(self._df)))

        all_properties = set(self._supported_position_properties + self._supported_sequence_properties + self._supported_mutation)

        diff = set(df['type']).difference(all_properties)
        if len(diff) > 0:
            self.log.warning("the following annotation types were not recognized: %s" % ", ".join(diff))

    def add_mutations(self, sequence, metadata=[]):

        if self._df is None:
            return

        tmp_df = self._df[ self._df['type'] == 'mutation' ]

        if 'genomic_mutations' in metadata:
            out_metadata = { 'genomic_mutations' : tmp_df['function'].str.split(',').tolist() }

        mutations = tmp_df['site'].tolist()

        unique_mutations = list(set(tmp_df['site']))

        self.log.info("unique mutations found in datafile: %s" % (", ".join(sorted(unique_mutations))))

        for m in unique_mutations:
            tokens = parse(self._mut_prot_parse, m)
            if tokens is None:
                self.log.error(f"Failed to parse mutation {m}; check that the format is correct (i.e. HGVS protein single amino acid substitution)")
                continue

            num = int(tokens['position'])
            wt  = three_to_one(str.upper(tokens['wt']))
            mut = three_to_one(str.upper(tokens['mut']))

            if wt == mut:
                self.log.info("synonymous mutation %s discarded" % m)
                continue

            try:
                site_seq_idx = sequence.seq2index(num)
            except:
                self.log.warning(f"mutation {m} is outside the protein sequence; it will be skipped")
                continue

            position = sequence.positions[site_seq_idx]

            if position.wt_residue_type != wt:
                self.log.warning("for mutation %s, residue %d is %s in wild-type sequence; it will be skipped" %(m, num, position.wt_residue_type))
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


    def add_position_properties(self, sequence, prop=_supported_position_properties):

        if self._df is None:
            return

        tmp_df = self._df[ self._df.apply(lambda x: x['type'] in self._supported_position_properties, axis=1) ]

        for idx, row in tmp_df.iterrows():

            try:
                this_position = sequence.positions[sequence.seq2index(int(row['site']))]
            except:
                self.log.error("position property refers to a position outside of the sequence")

            property_obj = position_properties_classes[row['type']](sources=[self],
                                                                    position=this_position)

            property_obj.metadata['function']  = row['function']
            property_obj.metadata['reference'] = row['reference']

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
                    self.log.error("line %d range in site column must be specified as 'from-to' (ex 10-15)" % (idx+1))
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
                                                                    name=row['name'])

            property_obj.metadata['function']  = row['function']
            property_obj.metadata['reference'] = row['reference']

            sequence.add_property(property_obj)
