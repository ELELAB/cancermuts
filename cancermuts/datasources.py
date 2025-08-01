# datasources.py - data sources handling for the cancermuts package
# (c) 2018 Matteo Tiberti <matteo.tiberti@gmail.com>
# (c) 2023 Katrine Meldgård <katrine@meldgaard.dk>
# (c) 2025 Pablo Sanchez-Izquierdo
# This file is part of cancermuts
# The function '_get_popmax_af' is taken and modified from the 'gnomad2csv' script
# which is part of the ELELAB/CSB-scripts repository
#
# cancermuts is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cancermuts is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with cancermuts.  If not, see <http://www.gnu.org/licenses/>.


"""
datasources classes --- :mod:`cancermuts.datasources`
================================================================
Classes to interrogate data sources and annotate various data

"""

import time
import requests as rq
from bioservices.uniprot import UniProt as bsUniProt
from Bio.PDB.Polypeptide import three_to_index, index_to_one
from Bio import SeqIO
import numpy as np
import pandas as pd
from .core import Sequence, Mutation
from .properties import *
from .metadata import *
from .log import *
from .exceptions import UnexpectedIsoformError
from io import StringIO
import json
import re
import sys
import os
from biothings_client import get_client
import pyliftover
import itertools
from parse import parse
from multiprocessing.dummy import Pool as threadPool
from future.utils import iteritems
from bravado.client import SwaggerClient
from bravado.requests_client import RequestsClient
from requests.adapters import HTTPAdapter
import gget
from io import StringIO


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

    def get_sequence(self, gene_id, upid=None, upac=None, isoform=None):
        if upac is not None:
            this_upac = upac
            self.log.info(f"The user-provided UniProt AC ({upac}) will be used")
        elif upid is not None:
            self.log.info("UniProt AC will be mapped from UniProt ID")
            this_upac = self._get_aliases(upid, ['UniProtKB_primaryAccession'])['UniProtKB_primaryAccession']
        else:
            self.log.info("retrieving UniProt ID for human gene %s" % gene_id)
            try:
                uniprot_table = pd.read_csv(StringIO(self._uniprot_service.search(f"(gene:{gene_id}) AND (organism_id:9606)")), sep='\t')
                upids = uniprot_table['Entry Name'].to_list()
            except:
                self.log.error("Failed to retrieve list of Uniprot IDs")
                return None

            this_upid = upids[0]

            if len(upids) > 1:
                self.log.warning("the following UniProt entries were found for gene %s: %s; will use %s" %(gene_id, ', '.join(upids), this_upid))
            else:
                self.log.info("will use Uniprot ID %s" % this_upid)

            this_upac = self._get_aliases(this_upid, ['UniProtKB_primaryAccession'])['UniProtKB_primaryAccession']

        if upid is None:
            this_upid = self._get_aliases(this_upac, ['UniProtKB_uniProtkbId'])['UniProtKB_uniProtkbId']
        else:
            this_upid = upid

        if isoform is not None:
            self.log.info(f"Isoform requested: {isoform}")

            url = f"https://rest.uniprot.org/uniprotkb/{this_upac}.json"
            response = rq.get(url)
            if not response.ok:
                raise ValueError(f"Failed to fetch UniProt JSON entry for {this_upac}")
            try:
                data = response.json()
            except Exception as e:
                raise ValueError(f"Invalid JSON response from UniProt for {this_upac}: {str(e)}")

            if "comments" not in data:
                raise ValueError(f"No 'comments' section found in UniProt entry for {this_upac}")

            alt_prods = [x for x in data["comments"] if x["commentType"] == "ALTERNATIVE PRODUCTS"]
            if not alt_prods:
                raise ValueError(f"No alternative products found in UniProt entry for {this_upac}. Cannot resolve isoform '{isoform}'.")

            is_canonical = False
            isoform_found = False
            try:
                for prod in alt_prods:
                    for iso in prod["isoforms"]:
                        ids = iso["isoformIds"]
                        status = iso["isoformSequenceStatus"]
                        if isoform in ids:
                            isoform_found = True
                            is_canonical = (status == "Displayed")
                            break
                    if isoform_found:
                        break
            except KeyError as e:
                raise ValueError(f"Missing expected field '{e.args[0]}' in ALTERNATIVE PRODUCTS for UniProt entry {this_upac}")

            if not isoform_found:
                raise ValueError(f"Isoform {isoform} not listed in UniProt entry {this_upac}")

            fasta_id = isoform
        else:
            fasta_id = this_upac
            is_canonical = True

        this_entrez = self._get_aliases(this_upac, ['GeneID'])

        if this_entrez is not None:
            this_entrez = this_entrez['GeneID']

            aliases = {'uniprot'     : this_upid,
                       'entrez'      : this_entrez,
                       'uniprot_acc' : this_upac }
        else:
            aliases = {'uniprot'     : this_upid,
                       'uniprot_acc' : this_upac }

        self.log.info("final aliases: %s" % aliases)

        self.log.info("retrieving sequence for UniProt sequence for Uniprot ID %s, Uniprot AC %s, gene %s" % (this_upid, fasta_id, gene_id))

        try:
            sequence = self._get_fasta_sequence(fasta_id)
        except:
            self.log.error("failed retrieving sequence for Uniprot ID %s" % fasta_id)
            return None

        return Sequence(gene_id=gene_id, uniprot_ac=this_upac, sequence=sequence, source=self, isoform=isoform, is_canonical=is_canonical, aliases=aliases)

    def _get_fasta_sequence(self, upid):
        """This function downloads the sequence of the specified UniProt entry
        by means of the UniProt API. The downloaded sequence corresponds to the
        sequence of the main isoform.

        Parameters
        ----------
        upid : :obj:`str`
            a valid UniProt identifier (e.g. Uniprot ID or Uniprot AC)

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
            self.log.info("fetching alias, fr=%s, to=%s, query=%s" % (fr, t_query, gene_id))
            responses = self._uniprot_service.mapping( fr = fr,
                                                       to = t_query,
                                                       query = gene_id )['results']

            if len(responses) == 0:
                self.log.warning(f"No {t} found for {gene_id}")
                return None

            if len(responses) > 1:
                self.log.warning(f"Multiple {t} found for {gene_id}. Found {t}: {', '.join(response['to'] for response in responses)}. No {t} will be assigned.")
                out[t] = None
                continue

            results = responses[0]

            if t_keyword is not None:
                self.log.info('using extracted keyword %s to parse results' % t_keyword)
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
        _cBioPortal_supported_metadata = ['cancer_type', 'cancer_study', 'genomic_coordinates', 'genomic_mutations']

        if not sequence.is_canonical:
            raise UnexpectedIsoformError("cBioPortal mutation annotation only supports canonical isoforms. Please use a Sequence object for a canonical isoform")

        if 'entrez' not in sequence.aliases.keys() or sequence.aliases['entrez'] is None:
            self.log.error('Entrez ID alias not available in sequence object')
            raise TypeError('Entrez ID alias not available in sequence object')

        for md in metadata:
            if md not in _cBioPortal_supported_metadata:
                self.log.error(f'{md} is not a valid metadata. Supported metadata are: {_cBioPortal_supported_metadata}')
                raise ValueError(f'{md} is not a valid metadata. Supported metadata are: {_cBioPortal_supported_metadata}')

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

            #Finds the mutation's indices in the complete list of mutations. Used for adding metadata
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
    def __init__(self, targeted_database_file, screen_mutant_database_file, classification_database_file, transcript_database_file, database_encoding=None, lazy_load_db=True):
        description = "COSMIC Database"
        super(COSMIC, self).__init__(name='COSMIC', version='v87', description=description)

        self._mut_regexp = 'p\.[A-Z][0-9]+[A-Z]$'
        self._mut_prog = re.compile(self._mut_regexp)
        self._mut_snv_regexp = '^[0-9]+:g\.[0-9+][ACTG]>[ACTG]'
        self._mut_snv_prog = re.compile(self._mut_snv_regexp)
        self._cosmic_phenotype_id_kwd = ['COSMIC_PHENOTYPE_ID']
        self._site_kwd = ['PRIMARY_SITE', 'SITE_SUBTYPE_1', 'SITE_SUBTYPE_2', 'SITE_SUBTYPE_3']
        self._histology_kwd = ['PRIMARY_HISTOLOGY', 'HISTOLOGY_SUBTYPE_1', 'HISTOLOGY_SUBTYPE_2', 'HISTOLOGY_SUBTYPE_3']

        self._use_cols_database_files = ['GENE_SYMBOL',
                                'TRANSCRIPT_ACCESSION',
                                'COSMIC_PHENOTYPE_ID',
                                 'MUTATION_AA',
                                 'CHROMOSOME',
                                 'GENOME_START',
                                 'GENOME_STOP',
                                 'MUTATION_CDS',
                                 'STRAND',
                                 'HGVSG']

        self._use_cols_classification_files = self._cosmic_phenotype_id_kwd + self._site_kwd + self._histology_kwd

        self._use_cols_transcript_files = ['TRANSCRIPT_ACCESSION','IS_CANONICAL']


        database_files = [targeted_database_file,screen_mutant_database_file,
                          classification_database_file,transcript_database_file]
        for file in database_files:
            if not isinstance(file, str):
                self.log.error('COSMIC database file must be a string.')
                raise TypeError('COSMIC database file  must be a string.')

        self._targeted_database_file = targeted_database_file
        self._screen_mutant_database_file = screen_mutant_database_file
        self._classification_database_file = classification_database_file
        self._transcript_database_file = transcript_database_file

        if database_encoding is None or isinstance(database_encoding, str):
            self._encoding = database_encoding
        else:
            self.log.errror('encoding for COSMIC database files must be None, or a single string that applies to all files')
            raise TypeError('encoding for COSMIC database files must be None, or a single string that applies to all files')

        for file in database_files:
            try:
                file = open(file, 'r')
                file.close()
            except Exception as e:
                self.log.error(f'Error in reading database file {e}')

        if not lazy_load_db:
            self._df = self._load_db_files(self._targeted_database_file, self._screen_mutant_database_file)
        else:
            self._df = None

    def _load_db_files(self, targeted_db_file, screenmut_db_file):

        targeted_screenmut__db_files = [targeted_db_file, screenmut_db_file]

        targeted_screenmut_dataframes = []
        for fi, file in enumerate(targeted_screenmut__db_files):
            self.log.info(f"Parsing database file {fi+1}...")
            try:
                targeted_screenmut_dataframes.append(pd.read_csv(file, sep='\t', dtype='str', na_values='NS', usecols=self._use_cols_database_files, encoding=self._encoding))
            except:
                self.log.error(f"Couldn't parse database file {fi+1}")
                raise TypeError(f"Couldn't parse database file {fi+1}")

        tmp_targeted_screenmut_df = pd.concat(targeted_screenmut_dataframes, ignore_index=True, sort=False)

        self.log.info(f"Parsing database file {len(targeted_screenmut__db_files)+1}...")
        try:
            classification_df = pd.read_csv(self._classification_database_file, sep='\t', dtype='str', na_values='NS', usecols=self._use_cols_classification_files, encoding=self._encoding)
        except ValueError:
            self.log.error(f"Couldn't parse database file {len(targeted_screenmut__db_files)+1}")
            raise TypeError(f"Couldn't parse database file {len(targeted_screenmut__db_files)+1}")

        self.log.info(f"Parsing database file {len(targeted_screenmut__db_files)+2}...")
        try:
            transcript_df = pd.read_csv(self._transcript_database_file, sep='\t', dtype='str', na_values='NS', usecols=self._use_cols_transcript_files, encoding=self._encoding)
        except ValueError:
            self.log.error(f"Couldn't parse database file {len(targeted_screenmut__db_files)+2}")
            raise TypeError(f"Couldn't parse database file {len(targeted_screenmut__db_files)+2}")


        self.log.info("Merging database files into a dataframe...")
        try:
            tmp_targeted_screenmut_classification_df = tmp_targeted_screenmut_df.merge(classification_df, on = 'COSMIC_PHENOTYPE_ID', sort=False)
        except KeyError:
            self.log.error("Couldn't merge database files due to missing or incorrectly named join columns")
            raise TypeError("Couldn't merge database files due to missing or incorrectly named join columns")

        try:
            df = tmp_targeted_screenmut_classification_df.merge(transcript_df, on = 'TRANSCRIPT_ACCESSION',sort=False)
        except KeyError:
            self.log.error("Couldn't merge database files due to missing or incorrectly named join columns")
            raise TypeError("Couldn't merge database files due to missing or incorrectly named join columns")

        return df

    def _parse_db_files(self, gene_id, genome_assembly_version = 'GRCh38',
                       cancer_types=None,
                       cancer_histology_subtype_1=None,
                       cancer_histology_subtype_2=None,
                       cancer_histology_subtype_3=None,
                       cancer_sites=None,
                       cancer_site_subtype_1=None,
                       cancer_site_subtype_2=None,
                       cancer_site_subtype_3=None,
                       metadata=[]):

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


        if self._df is not None:
            df = self._df[ (self._df['GENE_SYMBOL'] == gene_id) & (self._df['IS_CANONICAL'] == 'y') ]
        else:
            filtered_lines_t = []
            filtered_lines_s = []
            with open(self._targeted_database_file, "r", encoding=self._encoding) as t, open(self._screen_mutant_database_file, "r", encoding=self._encoding) as s:

                filtered_lines_t.append(t.readline())
                filtered_lines_s.append(s.readline())
                filtered_lines_t += [ line for line in t if gene_id in line.split() ]
                filtered_lines_s += [ line for line in s if gene_id in line.split() ]

            if len(filtered_lines_s) == 1 and len(filtered_lines_t) == 1:
                raise ValueError(f"The given gene_id {gene_id} is not present in the database files")

            filtered_targeted_database_file = StringIO("".join(filtered_lines_t))
            filtered_screenmut_database_file = StringIO("".join(filtered_lines_s))

            df = self._load_db_files(filtered_targeted_database_file, filtered_screenmut_database_file)

            df = df[ (df['IS_CANONICAL'] == 'y') ]

        if cancer_types is not None:
            df = df[ df['PRIMARY_HISTOLOGY'].isin(cancer_types) ]
        if cancer_histology_subtype_1 is not None:
            df = df[ df['HISTOLOGY_SUBTYPE_1'].isin(cancer_histology_subtype_1) ]
        if cancer_histology_subtype_2 is not None:
            df = df[ df['HISTOLOGY_SUBTYPE_2'].isin(cancer_histology_subtype_2) ]
        if cancer_histology_subtype_3 is not None:
            df = df[ df['HISTOLOGY_SUBTYPE_3'].isin(cancer_histology_subtype_3) ]

        if cancer_sites is not None:
            df = df[ df['PRIMARY_SITE'].isin(cancer_sites) ]
        if cancer_site_subtype_1 is not None:
            df = df[ df['SITE_SUBTYPE_1'].isin(cancer_site_subtype_1) ]
        if cancer_site_subtype_2 is not None:
            df = df[ df['SITE_SUBTYPE_2'].isin(cancer_site_subtype_2) ]
        if cancer_site_subtype_3 is not None:
            df = df[ df['SITE_SUBTYPE_3'].isin(cancer_site_subtype_3) ]

        df = df[ df['MUTATION_AA'].notna() ]

        df = df[ df.apply(lambda x: bool(self._mut_prog.match(x['MUTATION_AA'])), axis=1) ]

        for r in df.iterrows():
            r = r[1]
            mutations.append(r['MUTATION_AA'])

            if do_cancer_type:
                out_metadata['cancer_type'].append([r['PRIMARY_HISTOLOGY']])

            if do_genomic_coordinates or do_genomic_mutations:
                gd = []
                if not isinstance(genome_assembly_version, str):
                    raise TypeError(f"Incorrect format for genome assembly version")
                grch = str(genome_assembly_version)
                if grch == 'GRCh38'or grch == 'hg38':
                    gd.append('hg38')
                elif grch == 'GRCh37' or grch == 'hg19':
                    gd.append('hg19')
                else:
                    raise ValueError(f"Unsupported genome assembly version {grch}")


                gd.append(r['CHROMOSOME'])
                gd.append(r['GENOME_START'])
                gd.append(r['GENOME_STOP'])
                gd.append(r['MUTATION_CDS'][-3])

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



    def add_mutations(self, sequence, genome_assembly_version='GRCh38',
                    cancer_types=None,
                    cancer_histology_subtype_1=None,
                    cancer_histology_subtype_2=None,
                    cancer_histology_subtype_3=None,
                    cancer_sites=None,
                    cancer_site_subtype_1=None,
                    cancer_site_subtype_2=None,
                    cancer_site_subtype_3=None,
                    use_alias=None, metadata=[]):
        _cosmic_supported_metadata = ['cancer_type', 'genomic_coordinates', 'genomic_mutations', 'cancer_site', 'cancer_histology']

        for md in metadata:
            if md not in _cosmic_supported_metadata:
                self.log.error(f'{md} is not a valid metadata. Supported metadata are: {_cosmic_supported_metadata}')
                raise ValueError(f'{md} is not a valid metadata. Supported metadata are: {_cosmic_supported_metadata}')

        if cancer_types is None:
            self.log.info("no cancer type specified; will use all of them")
        if use_alias is not None:
            gene_id = sequence.aliases[use_alias]
            self.log.info("using alias %s as gene name" % sequence.aliases[use_alias])
        else:
            gene_id = sequence.gene_id

        raw_mutations, out_metadata = self._parse_db_files(gene_id, genome_assembly_version = genome_assembly_version,
                                                            cancer_types=cancer_types,
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

        if not sequence.is_canonical:
            raise UnexpectedIsoformError(
                "PhosphoSite annotation only supports canonical isoforms. Please use a Sequence object for a canonical isoform")

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
        description = "MyVariant.Info Database via BioThings client"
        super(MyVariant, self).__init__(name='MyVariant', version='1,0', description=description)

        self._biothings = get_client('variant')
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

        mutation.metadata['revel_score'] = []

        gcs = mutation.metadata.get('genomic_coordinates', [])
        gms = mutation.metadata.get('genomic_mutations', [])

        if (gcs is None or gcs == []) and (gms is None or gms == []):
            self.log.warning(f"No genomic coordinates or mutations for {mutation}; skipping.")
            return False

        if gcs is None or gcs == []:
            self.log.warning(f"No genomic coordinates available for revel score, {mutation}")
            gcs = [None] * len(gms)

        if gms is None or gcs == []:
            self.log.warning(f"No genomic mutations available for revel score, {mutation}")
            gcs = [None] * len(gms)


        gcs_gms = list(zip(gcs, gms))

        for gc,gm in gcs_gms:
            self.log.debug(f"Pulling revel score for {mutation}, {gc}, {gm}")
            if gm is not None:
                self.log.debug(f"Genomic mutation will be used to retrieve revel score for mutation {mutation}")
                revel_scores = self._get_revel_from_gm(mutation, gm)
            elif gc is not None:
                self.log.debug(f"Genomic coordinates will be used to retrieve revel score for mutation {mutation}")
                revel_scores = self._get_revel_from_gc(mutation, gc)

            self.log.debug(f"Final revel score for {gc} {gm} {revel_scores}")

            if revel_scores is not None:
                mutation.metadata['revel_score'].extend(revel_scores)

        self.log.debug(f"Final revel metadata for {mutation} {mutation.metadata['revel_score']}")

    def _validate_revel_hit(self, mutation, hit):

        try:
            aa = hit['dbnsfp']['aa']
        except:
            self.log.error("No residue information was found; it will be skipped.")
            return False

        if not isinstance(aa, list):
            aa = [aa]

        aa_short = []
        self.log.debug(f"{len(aa)} residue definitions will be tested")
        for idx_this_aa, this_aa in enumerate(aa):
            try:
                this_hit_pos = this_aa['pos']
                this_hit_ref = this_aa['ref']
                this_hit_alt = this_aa['alt']
            except:
                self.log.info(f"No residue information was found in variant definition {idx_this_aa}")
                continue

            if not mutation.sequence_position.wt_residue_type == this_hit_ref:
                self.log.info(f"Reference residue in revel does not correspond in variant definition {idx_this_aa}")
                continue

            if not mutation.sequence_position.sequence_position == this_hit_pos:
                try:
                    if not mutation.sequence_position.sequence_position in [ int(a) for a in this_hit_pos ]:
                        raise TypeError
                except:
                    self.log.info(f"Sequence position in revel does not correspond for in variant definition {idx_this_aa}")
                    continue

            if mutation.mutated_residue_type != this_hit_alt:
                self.log.info(f"Protein mutation in revel does not correspond for this in variant definition {idx_this_aa}")
                continue

            aa_short.append(idx_this_aa)

            if len(aa_short) == 1:
                self.log.info(f"A single valid residue entry was found for revel score in mutation {mutation}")
                return True
            elif len(aa_short) > 1:
                self.log.warning(f"More than one valid residue entry found for revel score in mutation {mutation}")
                return True
            elif len(aa_short) == 0:
                self.log.warning(f"No valid residue entry found for revel score in mutation {mutation}")
            return False

    def _convert_hg38_to_hg19(self, gc):
        if gc.genome_build == 'hg38':
            self.log.info("%s has genomic data in hg38 assembly - will be converted to hg19" % gc)
            converted_coords = self._lo.convert_coordinate('chr%s' % gc.chr, int(gc.get_coord()))

            if len(converted_coords) != 1:
                self.log.error("Could not convert genomic coordinates for %s (liftOver returned %d coords); it will be skipped" % (gc, len(converted_coords)))
                return None

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
            revel_scores = self._revel_cache_gm[query_str]
            self.log.info(f"Revel score for {query_str} retrieved from cache")
            return [ DbnsfpRevel(source=self, score=r) for r in revel_scores ]
        except KeyError:
            pass

        self.log.debug(f"Querying myvariant with string {query_str}.")
        hit = self._biothings.getvariant(query_str)

        if hit is None:
            self.log.error(f"Variant {query_str} was not found in MyVariant!")
            return None

        if not self._validate_revel_hit(mutation, hit):
            return None

        try:
            revel_scores = hit['dbnsfp']['revel']['score']
            self.log.debug(f"Downloaded revel score: {revel_scores}")
        except KeyError:
            self.log.warning("No revel score found for mutation {mutation} in this hit; it will be skipped")
            return None

        try:
            revel_scores = [float(revel_scores)]
        except TypeError:
            pass

        revel_scores = sorted(list(set(revel_scores)))
        self._revel_cache_gm[query_str] = revel_scores

        return [ DbnsfpRevel(source=self, score=r) for r in revel_scores ]

    def _get_revel_from_gc(self, mutation, gc):

        found_scores = []

        if gc is None:
            return None

        if gc.coord_start != gc.coord_end:
            self.log.warning(f"Mutation {mutation} has more than one nucleotide change; it will be skipped")
            return None

        converted_coords = self._convert_hg38_to_hg19(gc)
        if converted_coords is None:
            return None

        query_str = 'chr%s:%d' % converted_coords
        query = self._biothings.query(query_str)
        hits = query['hits']
        if len(hits) < 1:
            self.log.warning(f"No DBnsfp hits for mutation {mutation}; it will be skipped")
            return None
        else:
            self.log.info(f"{len(hits)} DBnsfp hits for mutation {mutation}")
        revel_scores = []

        for hit in hits:

            revel_score = None

            if not self._validate_revel_hit(mutation, hit):
                continue

            try:
                revel_score = hit['dbnsfp']['revel']['score']
            except KeyError:
                self.log.warning(f"No revel score found for mutation {mutation} in this hit; it will be skipped")
                continue

            try:
                revel_score = [float(revel_score)]
            except TypeError:
                pass

            revel_score = sorted(list(set(revel_score)))

            revel_scores.append(revel_score)

        if len(revel_scores) > 1: # more than one hit had usable scores
            self.log.warning("More than one hit had usable scores; it will be skipped")
            return None
        elif len(revel_scores) == 0:
            self.log.warning(f"No revel scores found using genomic coordinates for mutation {mutation}")
            return None
        else:
            revel_scores = revel_scores[0]
            self.log.debug(f"Downloaded revel score: {revel_scores}")
            return [ DbnsfpRevel(source=self, score=r) for r in revel_scores ]


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

    def add_sequence_properties(self, sequence, exclude_elm_classes=r'{.*}', use_alias='uniprot_acc', elm_wait=180):
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
                                                                         name=self._elm_classes[d[0]][0],
                                                                         id=d[0])

            property_obj.metadata['function'] = [self._elm_classes[d[0]][0]]
            property_obj.metadata['ref']      = self.description
            sequence.add_property(property_obj)

class ggetELMPredictions(StaticSource, object):
    @logger_init
    def __init__(self):
        """
        Data source for ELM which uses the gget Python package, rather than
        interrogating the ELM webserver
        """

        description = "ELM Prediction with gget"
        super(ggetELMPredictions, self).__init__(name='ggetELM', version='1.0', description=description)

    def _get_prediction(self, sequence):
        """
        Gets predicted SLIMs using regexp mode only, using the gget Python
        package

        Parameters
        ----------
        sequence : :obj:`str`
            Protein sequence, as a single string

        Returns
        ----------
        slims : :obj:`pandas.DataFrame`
            data frame containing SLIM definitions in form of:
              - ELM identifier
              - Name
              - Description
              - Start position
              - End position
        """

        try:
            ortho_slims, regex_slims = gget.elm(sequence, uniprot=False)
        except FileNotFoundError:
            gget.setup('elm')
            ortho_slims, regex_slims = gget.elm(sequence, uniprot=False)

        return regex_slims[['ELMIdentifier',
                            'FunctionalSiteName',
                            'Description',
                            'motif_start_in_query',
                            'motif_end_in_query']].drop_duplicates()

    def add_sequence_properties(self, sequence, exclude_elm_classes=r'{.*}'):
        """
        Adds sequence properties to a sequence object

        Parameters
        ----------
        sequence : :obj:`cancermuts.core.Sequence`
            Sequence object with the protein to be annotated

        exclude_elm_classes : :obj:`str`
            Regular expression matching ELM classes to be excluded from the output
        """

        self.log.info("adding gget ELM predictions to sequence ...")

        data = self._get_prediction(sequence.sequence)

        for _, r in data.iterrows():

            if re.match(exclude_elm_classes, r['ELMIdentifier']):
                self.log.info("%s was filtered out as requested" % r['ELMIdentifier'])
                continue

            this_positions = []
            for p in range(r['motif_start_in_query'], r['motif_end_in_query']+1):
                this_positions.append(sequence.positions[sequence.seq2index(p)])

            property_obj = sequence_properties_classes['linear_motif']  (sources=[self],
                                                                         positions=this_positions,
                                                                         name=r['FunctionalSiteName'],
                                                                         id=r['ELMIdentifier'])

            property_obj.metadata['function'] = [r['Description']]
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

    #Supported metadata for the different gnomAD versions
    _v2_1=['gnomad_exome_allele_frequency', 'gnomad_genome_allele_frequency', 'gnomad_popmax_exome_allele_frequency', 'gnomad_popmax_genome_allele_frequency']
    _v3=['gnomad_genome_allele_frequency', 'gnomad_popmax_genome_allele_frequency']

    _version_metadata_compatability = { '2.1' : _v2_1,
                                        '2.1_controls' : _v2_1,
                                        '2.1_non-neuro' : _v2_1,
                                        '2.1_non-cancer' : _v2_1,
                                        '2.1_non-topmed' : _v2_1,
                                        '3' : _v3
                                        }

    #Supported data for gnomAD versions
    _exome_genome_support = ['2.1', '2.1_controls', '2.1_non-neuro', '2.1_non-cancer', '2.1_non-topmed']
    _genome_support = ['3']

    @logger_init
    def __init__(self, version='2.1'):

        self._gnomad_version = str(version)
        if self._gnomad_version not in self._versions.keys():
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
        gnomADPopmaxExomeAlleleFrequency.set_version_in_desc(self._version_str[version])
        gnomADPopmaxGenomeAlleleFrequency.set_version_in_desc(self._version_str[version])

    def add_metadata(self, sequence, md_type=['gnomad_exome_allele_frequency'], use_alias=None):
        for md in md_type:
            if md not in self._version_metadata_compatability[self._gnomad_version]:
                self.log.error(f"The  metadata type '{md}' is not compatible with the gnomAD version. Compatible metadata are {self._version_metadata_compatability[self._gnomad_version]}")
                raise ValueError(f"The  metadata type '{md}' is not compatible with the gnomAD version. Compatible metadata are {self._version_metadata_compatability[self._gnomad_version]}")

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
        self.log.info("getting exome allele frequency")
        self._get_metadata(mutation, gene_id, 'gnomad_exome_allele_frequency')

    def _get_genome_allele_freq(self, mutation, gene_id):
        self.log.info("getting genome allele frequency")
        self._get_metadata(mutation, gene_id, 'gnomad_genome_allele_frequency')

    def _get_popmax_exome_allele_freq(self, mutation, gene_id):
        """Adds popmax exome allele frequency to mutation object.

        Parameters
        ----------
        mutation : :obj:`cancermuts.core.Mutation`
            Object containing mutations and metadata associated with the sequence.
        gene_id : :obj:`str`
            ID of the gene (gene name) to which the sequence belongs
        """

        self.log.info("getting popmax exome allele frequency")
        self._get_metadata(mutation, gene_id, 'gnomad_popmax_exome_allele_frequency')

    def _get_popmax_genome_allele_freq(self, mutation, gene_id):
        """Adds popmax genome allele frequency to mutation object.

        Parameters
        ----------
        mutation : :obj:`cancermuts.core.Mutation`
            Object containing mutations and metadata associated with the sequence.
        gene_id : :obj:`str`
            ID of the gene (gene name) to which the sequence belongs
        """

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
                try:
                    v_str = variant.as_assembly(ref_assembly).get_value_str(fmt='gnomad')
                except TypeError as e:
                    self.log.error(str(e))
                    v_str = None
                    continue
            else:
                v_str = None
                continue

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
                if pd.isna(this_df[exac_key[md_type]].values[0]):
                    af = None
                    mutation.metadata[md_type] = []
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

        if self._gnomad_version in self._exome_genome_support:
            variants['total_ac'] = variants['exome_ac'] + variants['genome_ac']
            variants['total_an'] = variants['exome_an'] + variants['genome_an']

            variants['exome_af']  =  variants['exome_ac'] / variants['exome_an']
            variants['genome_af'] = variants['genome_ac'] / variants['genome_an']
            variants['total_af'] = variants['total_ac'] / variants['total_an']

        elif self._gnomad_version in self._genome_support:
            variants['total_ac'] = pd.NA
            variants['total_an'] = pd.NA

            variants['exome_af']  =  pd.NA
            variants['genome_af'] = variants['genome_ac'] / variants['genome_an']
            variants['total_af'] = pd.NA

        variants[['popmax_exome_ac',  'popmax_exome_an',  'popmax_exome_af',
                'popmax_genome_ac', 'popmax_genome_an', 'popmax_genome_af',
                'popmax_tot_ac',    'popmax_tot_an',    'popmax_tot_af']] = variants.apply(self._get_popmax_af, axis=1)

        return variants

    def _get_popmax_af(self, variants, pops=['afr', 'eas', 'nfe', 'amr', 'sas']):
        """Calculates the popmax allele frequency of genome and/or exome, depending on the data available.

        Parameters
        ----------
        variants : :obj:`DataFrame`
            Dataframe containing at least the following:
                Index:
                    RangeIndex
                Columns:
                    'exome.populations' and/or 'genome.populations' :
                        list of dicts in the format [{'id': , 'ac': , 'an': }]
        pops : :obj:`list`
            Optional. List containing any or all of the populations
            ('afr', 'eas', 'nfe', 'amr', 'sas').

        Returns
        -------
        `Series`
            Series containing the calculated popmax values
        """
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
                if pd.isna(popmax_exome['af']).all():
                    popmax_exome_ac = pd.NA
                    popmax_exome_an = pd.NA
                    popmax_exome_af = pd.NA
                else:
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
                if pd.isna(popmax_genome['af']).all():
                    popmax_genome_ac = pd.NA
                    popmax_genome_an = pd.NA
                    popmax_genome_af = pd.NA
                else:
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
                        popmax_tot_ac,    popmax_tot_an,    popmax_tot_af])

class MobiDB(DynamicSource):

    description = "MobiDB"

    @logger_init
    def __init__(self):

        super(MobiDB, self).__init__(name='MobiDB', version='0.1', description=self.description)

        self._requests_url = 'https://mobidb.org/api/download?format=json&acc='
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

        # Collect data from MobiDB response
        try:
            curated_disorder_data = data['curated-disorder-priority']
        except KeyError:
            self.log.error("No curated disorder data was found.")
            curated_disorder_data = None

        try:
            derived_disorder_data = data['derived-disorder-priority']
        except KeyError:
            self.log.error("No derived disorder data was found.")
            derived_disorder_data = None

        try:
            homology_disorder_data = data['homology-disorder-priority']
        except KeyError:
            self.log.error("No homology based disorder data was found.")
            homology_disorder_data = None

        try:
            predicted_disorder_data = data['prediction-disorder-priority']
        except KeyError:
            self.log.error("No predicted disorder data was found.")
            predicted_disorder_data = None

        # Initialisation for annotation
        # Disorder types should be in order of priority, with the highest priority first
        disorder_types=['curated-disorder-priority','derived-disorder-priority','homology-disorder-priority','prediction-disorder-priority']
        evidence_types = {'curated-disorder-priority': 'curated',
                        'derived-disorder-priority': 'derived',
                        'homology-disorder-priority': 'homology',
                        'prediction-disorder-priority': 'prediction'}

        disorder_data={'curated-disorder-priority': curated_disorder_data,
                        'derived-disorder-priority': derived_disorder_data,
                        'homology-disorder-priority': homology_disorder_data,
                        'prediction-disorder-priority': predicted_disorder_data}

        assignments = [None for p in sequence.positions]

        # If specified data is available, annotate sequence accordingly
        for dis_term in disorder_types:
            if disorder_data[dis_term] is not None:
                regions = disorder_data[dis_term]['regions']
                for r in regions:
                    for i in range(r[0]-1, r[1]): # r[1]: -1 because of the 0-offset, +1 because of the [) of range, total 0
                        if assignments[i] is None:
                            try:
                                assignments[i] = 'Disordered, '+evidence_types[dis_term]
                            except IndexError:
                                self.log.error("residue index %s not in sequence!" % i)
                                return None

        return assignments

    def _get_mobidb(self, sequence, use_alias='uniprot_acc'):
        if use_alias is not None:
            gene_id = sequence.aliases[use_alias]
            self.log.info("using alias %s as gene name" % sequence.aliases[use_alias])
        else:
            gene_id = sequence.gene_id

        # Download request from MobiDB
        try:
            url = ('').join([self._requests_url, gene_id])
            self.log.debug("fetching %s" % url)
            req = rq.get(url)
        except:
            self.log.error("Couldn't get data for %s" % gene_id)
            return None

        #Check if MobiDB responded correctly
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
    _metadata_cols = ['genomic_mutations']
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
        try:
            self._metadata_df = df[ self._metadata_cols ]
        except:
            self.log.warning("Metadata columns not found in csv file (these are: %s)." % ", ".join(self._metadata_cols))

        self.log.info("Parsed file:")
        self.log.info('\n{0}'.format(str(self._df)))

        all_properties = set(self._supported_position_properties + self._supported_sequence_properties + self._supported_mutation)

        diff = set(df['type']).difference(all_properties)
        if len(diff) > 0:
            self.log.warning("the following annotation types were not recognized: %s" % ", ".join(diff))

    def add_mutations(self, sequence, metadata=[]):

        if self._df is None:
            return

        #Prepares metadata for mutations
        do_genomic_mutations=False
        if 'genomic_mutations' in metadata:
            do_genomic_mutations=True
            try:
                gm_df = self._metadata_df['genomic_mutations'][self._df['type'] == 'mutation']
            except:
                self.log.error('genomic_mutations specified in metadata, but no genomic_mutations column found in file. Metadata annotation will be skipped.')
                metadata.remove('genomic_mutations')
                do_genomic_mutations=False

        if do_genomic_mutations:
            gm = list()
            for row in gm_df:
                if row != '':
                    gm.append(re.split('[,\s]+',row))
                else:
                    gm.append(None)
            out_metadata = {'genomic_mutations':gm}

        #mutation dataframe
        tmp_df = self._df[ self._df['type'] == 'mutation' ]


        mutations = tmp_df['site'].tolist()

        unique_mutations = list(set(tmp_df['site']))

        self.log.info("unique mutations found in datafile: %s" % (", ".join(sorted(unique_mutations))))

        for m in unique_mutations:
            #collects the position, wild type aa and mutated aa
            tokens = parse(self._mut_prot_parse, m)
            if tokens is None:
                self.log.error(f"Failed to parse mutation {m}; check that the format is correct (i.e. HGVS protein single amino acid substitution)")
                continue

            #checks if it is a synonymous mutation
            num = int(tokens['position'])
            wt  = index_to_one(three_to_index(str.upper(tokens['wt'])))
            mut = index_to_one(three_to_index(str.upper(tokens['mut'])))

            if wt == mut:
                self.log.info("synonymous mutation %s discarded" % m)
                continue

            #Checks if mutation is outside protein sequence
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

            if do_genomic_mutations:
                #If there are multiple genomic mutation metadata for one mutation, move the
                #additional entries to the end of out_metadata['genomic_mutations'] and
                #and add the index number for this entry to mutation_indices
                new_mut_ind = list()
                for mi in mutation_indices:
                    if out_metadata['genomic_mutations'][mi] is not None and len(out_metadata['genomic_mutations'][mi]) > 2:
                        num_md = len(out_metadata['genomic_mutations'][mi])
                        for i in range(2,num_md,2):
                            out_metadata['genomic_mutations'].append(out_metadata['genomic_mutations'][mi][i:i+2])
                            new_mut_ind.append(len(out_metadata['genomic_mutations'])-1)
                        out_metadata['genomic_mutations'][mi] = out_metadata['genomic_mutations'][mi][0:2]
                mutation_indices.extend(new_mut_ind)

            mutation_obj = Mutation(sequence.positions[site_seq_idx],
                                    mut,
                                    [self])

            #Adding mutation to mutation object
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
