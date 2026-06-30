# metadata.py - metadata handling for the cancermuts package
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
metadata classes --- :mod:`cancermuts.metadata`
================================================================
Classes to handle metadata

"""

from .log import logger_init
import pyliftover
from parse import parse
import re

lo_hg38_hg19 = pyliftover.LiftOver('hg38', 'hg19')
lo_hg19_hg38 = pyliftover.LiftOver('hg19', 'hg38')

class Metadata(object):
    def __init__(self, source):
        self.source = source

class CancerType(Metadata):

    description = "Cancer type"
    header = "cancer_type"

    def __init__(self, source, cancer_type):
        super(CancerType, self).__init__(source)
        self.cancer_type = cancer_type


    def get_value(self):
        return self.cancer_type

    def get_value_str(self):
        return self.get_value()

    def __repr__(self):
        return "<CancerType %s from %s>" % (self.cancer_type, self.source.name)

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        return self.source == other.source and self.cancer_type == other.cancer_type

    def __hash__(self):
        return hash((self.source, self.cancer_type))

class CancerStudy(Metadata):

    description = "Cancer study"
    header = "cancer_study"

    def __init__(self, source, study_id):
        super(CancerStudy, self).__init__(source)
        self.study_id = study_id

    def get_value(self):
        return self.study_id

    def get_value_str(self):
        return self.get_value()

    def __repr__(self):
        return "<CancerStudy %s from %s>" % (self.study_id, self.source.name)

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        return self.source == other.source and self.study_id == other.study_id

    def __hash__(self):
        return hash((self.source, self.study_id))

class GenomicCoordinates(Metadata):

    description = "Genomic coordinates"
    header = "genomic_coordinates"

    def __init__(self, source, genome_build, chromosome, coord_start, coord_end, ref):
        super(GenomicCoordinates, self).__init__(source)
        self.genome_build = genome_build

        if chromosome == '23':
            chromosome = 'X'
        if chromosome == '24':
            chromosome = 'Y'

        self.chr = chromosome
        self.coord_start = coord_start
        self.coord_end = coord_end
        self.ref = ref

    def get_value(self):
        return [self.genome_build, self.chr, self.coord_start, self.coord_end, self.ref]

    def get_value_str(self):
        return "%s,%s:%s-%s" % (self.genome_build, self.chr, self.coord_start, self.coord_end)

    def get_coord(self):
        return self.coord_start

    def __repr__(self):
        return "<GenomicCoordinates %s:%s-%s in %s from %s>" % (self.chr, self.coord_start, self.coord_end, self.genome_build, self.source.name)

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        return self.source == other.source and \
        self.chr == other.chr and \
        self.coord_start == other.coord_start and \
        self.coord_end == other.coord_end and \
        self.ref == other.ref and \
        self.genome_build == other.genome_build

    def __hash__(self):
        return hash((self.source, self.chr, self.coord_start, self.coord_end, self.ref, self.genome_build))

class GenomicMutation(Metadata):

    description = "Genomic mutation"
    header = "genomic_mutation"

    _mut_snv_prog = re.compile(r'^(?P<chr>[0-9XY]+):g\.(?P<start>[0-9]+)(?P<ref>[ACTG])>(?P<alt>[ACTG])$')
    _mut_delins_prog = re.compile(r'^(?P<chr>[0-9XY]+):g\.(?P<start>[0-9]+)(?:_(?P<end>[0-9]+))?delins(?P<alt>[ACTG]+)$')
    _mut_inv_prog = re.compile(r'^(?P<chr>[0-9XY]+):g\.(?P<start>[0-9]+)_(?P<end>[0-9]+)inv$')
    _mut_del_prog = re.compile(r'^(?P<chr>[0-9XY]+):g\.(?P<start>[0-9]+)(?:_(?P<end>[0-9]+))?del$')
    _mut_ins_prog = re.compile(r'^(?P<chr>[0-9XY]+):g\.(?P<start>[0-9]+)_(?P<end>[0-9]+)ins(?P<alt>[ACTG]+)$')
    _mut_dup_prog = re.compile(r'^(?P<chr>[0-9XY]+):g\.(?P<start>[0-9]+)(?:_(?P<end>[0-9]+))?dup$')
    _mut_repeat_prog = re.compile(r'^(?P<chr>[0-9XY]+):g\.(?P<start>[0-9]+)(?:_(?P<end>[0-9]+))?(?P<repeat_unit>[ACTG]+)\[(?P<repeat_count>[0-9]+)\]$')

    @staticmethod
    def _normalise_chr(chromosome):
        if chromosome == '23':
            return 'X'
        elif chromosome == '24':
            return 'Y'
        return chromosome

    @logger_init
    def __init__(self, source, genome_build, definition):
        super(GenomicMutation, self).__init__(source)

        self.genome_build = genome_build
        self.definition = definition

        self.chr = None
        self.start = None
        self.end = None
        self.ref = None
        self.alt = None
        self.mutation_type = None
        self._parse_definition(definition)

    def _set_common_fields(self, tokens, mutation_type):
        self.mutation_type = mutation_type
        self.chr = self._normalise_chr(tokens["chr"])
        self.start = int(tokens["start"])
        self.end = int(tokens["end"]) if tokens.get("end") is not None else self.start
        self.ref = tokens.get("ref")
        self.alt = tokens.get("alt")

    def _position_string(self):
        if self.start == self.end:
            return str(self.start)
        else:
            return f"{self.start}_{self.end}"

    def _parse_definition(self, definition):
        match = self._mut_snv_prog.match(definition)
        if match:
            self._set_common_fields(match.groupdict(), "snv")
            self.definition = f"{self.chr}:g.{self.start}{self.ref}>{self.alt}"
            return

        match = self._mut_delins_prog.match(definition)
        if match:
            self._set_common_fields(match.groupdict(), "delins")
            self.definition = f"{self.chr}:g.{self._position_string()}delins{self.alt}"
            return

        match = self._mut_inv_prog.match(definition)
        if match:
            self._set_common_fields(match.groupdict(), "inversion")
            self.definition = f"{self.chr}:g.{self._position_string()}inv"
            return

        match = self._mut_ins_prog.match(definition)
        if match:
            self._set_common_fields(match.groupdict(), "insertion")
            self.definition = f"{self.chr}:g.{self.start}_{self.end}ins{self.alt}"
            return

        match = self._mut_dup_prog.match(definition)
        if match:
            self._set_common_fields(match.groupdict(), "duplication")
            self.definition = f"{self.chr}:g.{self._position_string()}dup"
            return

        match = self._mut_del_prog.match(definition)
        if match:
            self._set_common_fields(match.groupdict(), "deletion")
            self.definition = f"{self.chr}:g.{self._position_string()}del"
            return

        match = self._mut_repeat_prog.match(definition)
        if match:
            tokens = match.groupdict()
            tokens["alt"] = f"{tokens['repeat_unit']}[{tokens['repeat_count']}]"
            self._set_common_fields(tokens, "repeat")
            self.definition = f"{self.chr}:g.{self._position_string()}{self.alt}"
            return

        self.log.info(f"unsupported genomic mutation format: {definition}")

    def get_value_str(self, fmt='csv'):
        if fmt == 'csv':
            return f"{self.genome_build},{self.definition}"
        if fmt == 'gnomad':
            if self.mutation_type == "snv":
                return f"{self.chr}-{self.start}-{self.ref}-{self.alt}"
            return None

    def get_coord(self):
        return self.start

    def _as_build(self, target_build):
        if self.genome_build == target_build:
            return self

        if self.genome_build == 'hg38' and target_build == 'hg19':
            lo = lo_hg38_hg19
        elif self.genome_build == 'hg19' and target_build == 'hg38':
            lo = lo_hg19_hg38
        else:
            raise TypeError

        converted_coords_start = lo.convert_coordinate('chr%s' % self.chr, int(self.start))
        converted_coords_end = lo.convert_coordinate('chr%s' % self.chr, int(self.end))

        if len(converted_coords_start) != 1:
            raise TypeError("Could not convert genomic coordinates (liftOver returned %d coords)" % len(converted_coords_start))
        if len(converted_coords_end) != 1:
            raise TypeError("Could not convert genomic coordinates (liftOver returned %d coords)" % len(converted_coords_end))

        converted_start = converted_coords_start[0][1]
        converted_end = converted_coords_end[0][1]

        if converted_start == converted_end:
            position = str(converted_start)
        else:
            position = f"{converted_start}_{converted_end}"

        if self.mutation_type == "snv":
            definition = f"{self.chr}:g.{converted_start}{self.ref}>{self.alt}"

        elif self.mutation_type == "delins":
            definition = f"{self.chr}:g.{position}delins{self.alt}"

        elif self.mutation_type == "deletion":
            definition = f"{self.chr}:g.{position}del"

        elif self.mutation_type == "insertion":
            if abs(converted_end - converted_start) != 1:
                raise TypeError("Could not convert insertion coordinates: converted positions are not adjacent")
            definition = f"{self.chr}:g.{converted_start}_{converted_end}ins{self.alt}"

        elif self.mutation_type == "duplication":
            definition = f"{self.chr}:g.{position}dup"

        elif self.mutation_type == "inversion":
            definition = f"{self.chr}:g.{position}inv"

        elif self.mutation_type == "repeat":
            definition = f"{self.chr}:g.{position}{self.alt}"

        else:
            raise TypeError

        return GenomicMutation(self.source, target_build, definition)

    def as_hg19(self):
        return self._as_build('hg19')

    def as_hg38(self):
        return self._as_build('hg38')

    def as_assembly(self, assembly):
        if assembly == 'hg19' or assembly == 'GRCh37':
            return self.as_hg19()
        elif assembly == 'hg38' or assembly == 'GRCh38':
            return self.as_hg38()
        else:
            raise TypeError

    def __repr__(self):
        return "<GenomicMutation %s from %s>" % (self.get_value_str(), self.source.name)

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):

        return self.source == other.source and \
        self.definition == other.definition

    def __hash__(self):
        return hash((self.source, self.definition))

class Revel(Metadata):

    description = "Revel score"
    header = "REVEL_score"

    def __init__(self, source, score):
        super(Revel, self).__init__(source)
        self.source = source
        self.score = score

    def get_value(self):
        return self.score

    def get_value_str(self):
        return "%.3f" % self.score

    def __repr__(self):
        return "<Revel, %.3f>" % self.score

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        return self.source == other.source and \
                self.score == other.score

    def __hash__(self):
        return hash((self.source, self.score))

class gnomADAlleleFrequency(Metadata):

    freq_type = 'generic'

    def __init__(self, source, frequency):
        super(gnomADAlleleFrequency, self).__init__(source)
        self.source = source
        self.frequency = frequency

    def get_value(self):
        return self.frequency

    def get_value_str(self):
        if self.frequency is None:
            return ''
        return "%.16f" % self.frequency

    def __eq__(self, other):
        return self.source == other.source and \
               self.frequency == other.frequency and \
               self.freq_type == other.freq_type

    def __hash__(self):
        return hash((self.source, self.frequency, self.freq_type))

    def __str__(self):
        return self.__repr__()

class gnomADExomeAlleleFrequency(gnomADAlleleFrequency):

    freq_type = "Exome"
    basic_description = "%s allele frequency" % freq_type
    description = basic_description
    header = "gnomad_exome_af"

    @classmethod
    def set_version_in_desc(cls, version_string):
        cls.description = "%s (%s)" % (cls.basic_description, version_string)

    def __init__(self, source, frequency):
        super(gnomADExomeAlleleFrequency, self).__init__(source, frequency)

    def __repr__(self):
        return "<gnomADExomeAlleleFrequency, %f>" % self.frequency

class gnomADGenomeAlleleFrequency(gnomADAlleleFrequency):

    freq_type = "Genome"
    basic_description = "%s allele frequency" % freq_type
    description = basic_description
    header = "gnomad_genome_af"

    @classmethod
    def set_version_in_desc(cls, version_string):
        cls.description = "%s (%s)" % (cls.basic_description, version_string)

    def __init__(self, source, frequency):
        super(gnomADGenomeAlleleFrequency, self).__init__(source, frequency)

    def __repr__(self):
        return "<gnomADGenomeAlleleFrequency, %f>" % self.frequency

class gnomADPopmaxExomeAlleleFrequency(gnomADAlleleFrequency):
    """This class represents gnomad popmax exome allele frequency

    Attributes
    ----------
    basic_description : :obj:`str`
        string containing description of the frequency type
    description : :obj:`str`
        string containing description of the frequency type including gnomad version
    freq_type : :obj:`str``
        string containing the type of allele frequency
    get_value : :obj:`cancermuts.metadata.gnomADAlleleFrequency.get_value`
        gets the frequency from self
    get_value_str : :obj:`cancermuts.metadata.gnomADAlleleFrequency.get_value_str`
        gets the frequency from self and converts it to string
    header : :obj:`str``
        string containing the name of the column header for the metatable
    """

    freq_type = "Popmax exome"
    basic_description = "%s allele frequency" % freq_type
    description = basic_description
    header = "gnomad_popmax_exome_af"

    @classmethod
    def set_version_in_desc(cls, version_string):
        """Makes it possible to use the class for different versions

        Parameters
        ----------
        version_string
        """
        cls.description = "%s (%s)" % (cls.basic_description, version_string)

    def __init__(self, source, frequency):
        """Constructor for the gnomAD popmax exome allele frequency

        Parameters
        ----------
        source: :obj:`cancermuts.datasources.gnomAD`
            Object containing the source of the frequency
        version_string: 'str'
            String containing the gnomAD version
        """
        super(gnomADPopmaxExomeAlleleFrequency, self).__init__(source, frequency)

    def __repr__(self):
        """Function which creates the string to add to the mutation object

        Returns
        -------
        'str'
            String containing allele frequency type and the corresponding frequency.
        """
        return "<gnomADPopmaxExomeAlleleFrequency, %f>" % self.frequency

class gnomADPopmaxGenomeAlleleFrequency(gnomADAlleleFrequency):
    """This class represents gnomad popmax genome allele frequency

    Attributes
    ----------
    basic_description : :obj:`str`
        string containing description of the frequency type
    description : :obj:`str`
        string containing description of the frequency type including gnomad version
    freq_type : :obj:`str`
        string containing the type of allele frequency
    get_value : :obj:`cancermuts.metadata.gnomADAlleleFrequency.get_value`
        gets the frequency from self
    get_value_str : :obj:`cancermuts.metadata.gnomADAlleleFrequency.get_value_str`
        gets the frequency from self and converts it to string
    header : :obj:`str`
        string containing the name of the column header for the metatable
    """

    freq_type = "Popmax genome"
    basic_description = "%s allele frequency" % freq_type
    description = basic_description
    header = "gnomad_popmax_genome_af"

    @classmethod
    def set_version_in_desc(cls, version_string):
        """Makes it possible to use the class for different versions

        Parameters
        ----------
        version_string
        """
        cls.description = "%s (%s)" % (cls.basic_description, version_string)

    def __init__(self, source, frequency):
        """Constructor for the gnomAD popmax exome allele frequency

        Parameters
        ----------
        source: :obj:`cancermuts.datasources.gnomAD`
            Object containing the source of the frequency
        version_string: 'str'
            String containing the gnomAD version
        """
        super(gnomADPopmaxGenomeAlleleFrequency, self).__init__(source, frequency)

    def __repr__(self):
        """Function which creates the string to add to the mutation object

        Returns
        -------
        'str'
            String containing allele frequency type and the corresponding frequency.
        """
        return "<gnomADPopmaxGenomeAlleleFrequency, %f>" % self.frequency

class CancerSite(Metadata):

    description = "Cancer site"
    header = "cancer_site"
    not_specified = 'NS'

    def __init__(self, source, *sites):
        super(CancerSite, self).__init__(source)
        self.source = source
        sites = [ i for i in sites if i != self.not_specified ]
        sites = [ i for i in sites if i == i ] # remove nans
        self.sites = sites

    def get_value(self):
        return self.sites

    def get_value_str(self):
        return ", ".join(self.sites)

    def __repr__(self):
        return "<CancerSite, %s>" % self.get_value_str()

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        return self.source == other.source and \
               self.sites == other.sites

    def __hash__(self):
        return hash((self.source, self.site))

class CancerHistology(Metadata):

    description = "Cancer histology"
    header = "cancer_histology"
    not_specified = 'NS'

    def __init__(self, source, *histology):
        super(CancerHistology, self).__init__(source)
        histology = [ i for i in histology if i != self.not_specified ]
        histology = [ i for i in histology if i == i ] # remove nans
        self.source = source
        self.histology = histology

    def get_value(self):
        return self.histology

    def get_value_str(self):
        return ", ".join(self.histology)

    def __repr__(self):
        return "<TumorHistology, %s>" % self.get_value_str()

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        return self.source == other.source and \
               self.histology == other.histology

    def __hash__(self):
        return hash((self.source, self.histology))

class ClinVarVariantID(Metadata):
    description = "ClinVar variant ID"
    header = "clinvar_variant_id"

    def __init__(self, source, variant_id):
        super().__init__(source)
        self.variant_id = str(variant_id)

    def get_value(self):
        return self.variant_id

    def get_value_str(self):
        return self.variant_id

    def __repr__(self):
        return f"<ClinVarVariantID {self.variant_id} from {self.source.name}>"

class ClinvarClassification(Metadata):

    xml_key = ""
    description = ""
    header = ""

    def __init__(self, source, classification):
        super().__init__(source)
        self.classification = classification

    def get_value(self):
        return self.classification

    def get_value_str(self):
        return self.classification
    def __repr__(self):
       return "%s(source=%r, classification=%r)" % (self.__class__.__name__, self.source, self.classification)

class ClinvarGermlineClassification(ClinvarClassification):
    description = "Clinvar Germline classification"
    header = "clinvar_germline_classification"
    xml_key = "GermlineClassification"

class ClinvarClinicalImpactClassification(ClinvarClassification):
    description = "Clinvar Clinical impact classification"
    header = "clinvar_clinical_impact_classification"
    xml_key = "SomaticClinicalImpact"

class ClinvarOncogenicityClassification(ClinvarClassification):
    description = "Clinvar Oncogenicity classification"
    header = "clinvar_oncogenicity_classification"
    xml_key = "OncogenicityClassification"

class ClinvarReviewStatus(Metadata):

    xml_key = ""
    description = ""
    header = ""

    stars_map = {
        'practice guideline': 4,
        'reviewed by expert panel': 3,
        'criteria provided, multiple submitters, no conflicts': 2,
        'criteria provided, multiple submitters': 'NA',
        'criteria provided, conflicting interpretations': 1,
        'criteria provided, conflicting classifications': 1,
        'criteria provided, single submitter': 1,
        'no assertion for the individual variant': 0,
        'no interpretation for the single variant': 0,
        'no assertion criteria provided': 0,
        'no assertion provided': 0,
        'no classification provided': 0,
        'no classifications from unflagged records': 0
    }

    def __init__(self, source, status):
        super().__init__(source)
        self.status = status
        self.stars = self.stars_map[status]

    def get_value(self):
        return self.stars

    def get_value_str(self):
        return self.stars

    def get_status_text(self):
        return self.status

    def get_stars(self):
        return self.stars

    def __repr__(self):
       return "%s(source=%r, status=%r, stars=%r)" % (self.__class__.__name__, self.source, self.status, self.stars)

class ClinvarClinicalImpactReviewStatus(ClinvarReviewStatus):
    description = "Clinvar Clinical impact review status"
    header = "clinvar_clinical_impact_review_status"
    xml_key = "SomaticClinicalImpact"

    stars_map = {
        'practice guideline': 4,
        'reviewed by expert panel': 3,
        'criteria provided, multiple submitters': 2,
        'criteria provided, single submitter': 1,
        'no assertion criteria provided': 0,
        'no classification provided': 0,
        'no classification for the individual variant': 0,
        'no classifications from unflagged records': 0
    }

class ClinvarGermlineReviewStatus(ClinvarReviewStatus):
    description = "Clinvar Germline review status"
    header = "clinvar_germline_review_status"
    xml_key = "GermlineClassification"

class ClinvarOncogenicityReviewStatus(ClinvarReviewStatus):
    description = "Clinvar Oncogenicity review status"
    header = "clinvar_oncogenicity_review_status"
    xml_key = "OncogenicityClassification"

class ClinvarCondition(Metadata):

    xml_key = ""
    description = ""
    header = ""

    def __init__(self, source, conditions):
        super().__init__(source)
        self.conditions = conditions if isinstance(conditions, list) else [conditions]

    def get_value(self):
        return self.conditions

    def get_value_str(self):
        return ";".join(self.conditions)

    def __repr__(self):
       return "%s(source=%r, conditions=%r)" % (self.__class__.__name__, self.source, self.conditions)

class ClinvarGermlineCondition(ClinvarCondition):
    description = "Clinvar Germline condition"
    header = "clinvar_germline_condition"
    xml_key = "GermlineClassification"

class ClinvarClinicalImpactCondition(ClinvarCondition):
    description = "Clinvar Clinical impact condition"
    header = "clinvar_clinical_impact_condition"
    xml_key = "SomaticClinicalImpact"

class ClinvarOncogenicityCondition(ClinvarCondition):
    description = "Clinvar Oncogenicity condition"
    header = "clinvar_oncogenicity_condition"
    xml_key = "OncogenicityClassification"

metadata_classes = {
                     'cancer_type'                 : CancerType,
                     'cancer_study'                : CancerStudy,
                     'genomic_coordinates'         : GenomicCoordinates,
                     'genomic_mutations'           : GenomicMutation,
                     'revel_score'                 : Revel,
                     'gnomad_genome_allele_frequency' : gnomADGenomeAlleleFrequency,
                     'gnomad_exome_allele_frequency' : gnomADExomeAlleleFrequency,
                     'gnomad_popmax_genome_allele_frequency' : gnomADPopmaxGenomeAlleleFrequency,
                     'gnomad_popmax_exome_allele_frequency' : gnomADPopmaxExomeAlleleFrequency,
                     'cancer_site'                 : CancerSite,
                     'cancer_histology'            : CancerHistology,
                     'clinvar_variant_id'          : ClinVarVariantID,
                     'clinvar_germline_condition'          : ClinvarGermlineCondition,
                     'clinvar_germline_review_status'      : ClinvarGermlineReviewStatus,
                     'clinvar_germline_classification'     : ClinvarGermlineClassification,
                     'clinvar_oncogenicity_condition'      : ClinvarOncogenicityCondition,
                     'clinvar_oncogenicity_review_status'  : ClinvarOncogenicityReviewStatus,
                     'clinvar_oncogenicity_classification' : ClinvarOncogenicityClassification,
                     'clinvar_clinical_impact_condition'      : ClinvarClinicalImpactCondition,
                     'clinvar_clinical_impact_review_status'  : ClinvarClinicalImpactReviewStatus,
                     'clinvar_clinical_impact_classification' : ClinvarClinicalImpactClassification
                    }
