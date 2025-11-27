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

    allowed_bases = set(['A', 'C', 'G', 'T'])

    _mut_snv_regexp = '^[0-9XY]+:g\.[0-9]+[ACTG]>[ACTG]$'
    _mut_insdel_regexp = '^[0-9XY]+:g\.[0-9]+_[0-9]+delins[ACTG]+$'
    _mut_inv_regexp = '^[0-9XY]+:g\.[0-9]+_[0-9]+inv$'
    _mut_snv_prog = re.compile(_mut_snv_regexp)
    _mut_insdel_prog = re.compile(_mut_insdel_regexp)
    _mut_inv_prog = re.compile(_mut_inv_regexp)
    _mut_snv_parse = '{chr}:g.{coord:d}{ref:l}>{alt:l}'
    _mut_insdel_parse = '{chr}:g.{coord_start:d}_{coord_end:d}delins{substitution}'
    _mut_inv_parse = '{chr}:g.{coord_start:d}_{coord_end:d}inv'

    @logger_init
    def __init__(self, source, genome_build, definition):
        super(GenomicMutation, self).__init__(source)

        self.genome_build = genome_build
        self.definition = definition

        if self._mut_snv_prog.match(definition):
            tokens = parse(self._mut_snv_parse, definition)

            if tokens['chr'] == '23':
                self.chr = 'X'
            elif tokens['chr'] == '24':
                self.chr = 'Y'
            else:
                self.chr = tokens['chr']

            self.coord = tokens['coord']
            self.ref = tokens['ref']
            self.alt = tokens['alt']
            self.is_snv = True
            self.is_insdel = False
            self.is_inversion = False

            self.definition=f"{self.chr}:g.{self.coord}{self.ref}>{self.alt}"

        elif self._mut_insdel_prog.match(definition):
            tokens = parse(self._mut_insdel_parse, definition)

            if tokens['chr'] == '23':
                self.chr = 'X'
            elif tokens['chr'] == '24':
                self.chr = 'Y'
            else:
                self.chr = tokens['chr']

            self.coord_start = tokens['coord_start']
            self.coord_end = tokens['coord_end']
            self.substitution = tokens['substitution']
            self.is_snv = False
            self.is_insdel = True
            self.is_inversion = False

            self.definition = f"{self.chr}:g.{self.coord_start}_{self.coord_end}delins{self.substitution}"

        elif self._mut_inv_prog.match(definition):
            tokens = parse(self._mut_inv_parse, definition)

            if tokens['chr'] == '23':
                self.chr = 'X'
            elif tokens['chr'] == '24':
                self.chr = 'Y'
            else:
                self.chr = tokens['chr']

            self.coord_start = tokens['coord_start']
            self.coord_end = tokens['coord_end']
            self.substitution = None
            self.is_snv = False
            self.is_insdel = False
            self.is_inversion = True

            self.definition = f"{self.chr}:g.{self.coord_start}_{self.coord_end}inv"

        else:
            self.log.info("doing other")
            self.chr = None
            self.coord = None
            self.ref = None
            self.alt = None
            self.is_snv = False
            self.is_insdel = False
            self.is_inversion = False

    def get_value_str(self, fmt='csv'):
        if fmt == 'csv':
            return f"{self.genome_build},{self.definition}"
        if fmt == 'gnomad':
            if self.is_snv:
                return f"{self.chr}-{self.coord}-{self.ref}-{self.alt}"
            elif self.is_insdel:
                return f"{self.chr}-{self.coord_start}-?-{self.substitution}"
            else:
                return None
        else:
            return None

    def get_coord(self):
        return self.coord

    def as_hg19(self):
        if self.genome_build == 'hg19':
            return self
        elif self.genome_build == 'hg38':
            if self.is_snv:
                converted_coords = lo_hg38_hg19.convert_coordinate('chr%s' % self.chr, int(self.coord))
                if len(converted_coords) != 1:
                    raise TypeError("Could not convert genomic coordinates (liftOver returned %d coords)" % len(converted_coords))
                return GenomicMutation(self.source, 'hg19', f"{self.chr}:g.{converted_coords[0][1]}{self.ref}>{self.alt}")
            elif self.is_insdel:
                converted_coords_start = lo_hg38_hg19.convert_coordinate('chr%s' % self.chr, int(self.coord_start))
                converted_coords_end   = lo_hg38_hg19.convert_coordinate('chr%s' % self.chr, int(self.coord_end))
                if len(converted_coords_start) != 1:
                    raise TypeError("Could not convert genomic coordinates (liftOver returned %d coords)" % len(converted_coords_start))
                if len(converted_coords_end) != 1:
                    raise TypeError("Could not convert genomic coordinates (liftOver returned %d coords)" % len(converted_coords_end))

                return GenomicMutation(self.source, 'hg19', f"{self.chr}:g.{converted_coords_start[0][1]}_{converted_coords_end[0][1]}delins{self.substitution}")
            else:
                raise TypeError

    def as_hg38(self):
        if self.genome_build == 'hg38':
            return self
        elif self.genome_build == 'hg19':
            if self.is_snv:
                converted_coords = lo_hg19_hg38.convert_coordinate('chr%s' % self.chr, int(self.coord))
                if len(converted_coords) != 1:
                    raise TypeError("Could not convert genomic coordinates (liftOver returned %d coords)" % len(converted_coords))
                return GenomicMutation(self.source, 'hg38', f"{self.chr}:g.{converted_coords[0][1]}{self.ref}>{self.alt}")
            elif self.is_insdel:
                converted_coords_start = lo_hg19_hg38.convert_coordinate('chr%s' % self.chr, int(self.coord_start))
                converted_coords_end   = lo_hg19_hg38.convert_coordinate('chr%s' % self.chr, int(self.coord_end))
                if len(converted_coords_start) != 1:
                    raise TypeError("Could not convert genomic coordinates (liftOver returned %d coords)" % len(converted_coords_start))
                if len(converted_coords_end) != 1:
                    raise TypeError("Could not convert genomic coordinates (liftOver returned %d coords)" % len(converted_coords_end))
                return GenomicMutation(self.source, 'hg38', f"{self.chr}:g.{converted_coords_start[0][1]}{converted_coords_end[0][1]}delins{self.substitution}")
            else:
                raise TypeError

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

class ClinVarCondition(Metadata):
    description = "ClinVar condition"
    header = "clinvar_condition"

    def __init__(self, source, data):
        super().__init__(source)
        conds = data["GermlineClassification"]
        self.conditions = conds if isinstance(conds, list) else [conds]

    def get_value(self):
        return self.conditions

    def get_value_str(self):
        return ";".join(self.conditions)

class ClinVarReviewStatus(Metadata):
    description = "ClinVar review status"
    header = "clinvar_review_status"

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
        'no classification provided': 0
    }

    def __init__(self, source, data):
        super().__init__(source)
        self.status = data["GermlineClassification"]
        self.stars = self.stars_map[self.status]

    def get_value(self):
        return self.stars

    def get_value_str(self):
        return self.stars

    def get_status_text(self):
        return self.status

    def get_stars(self):
        return self.stars

    def __repr__(self):
        return f"<ClinVarReviewStatus {self.stars} (from: {self.status}) from {self.source.name}>"

class ClinVarClassification(Metadata):
    description = "ClinVar classification"
    header = "clinvar_classification"

    def __init__(self, source, data):
        super().__init__(source)
        self.classification = data["GermlineClassification"]

    def get_value(self):
        return self.classification

    def get_value_str(self):
        return self.classification

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
                     'clinvar_variant_id'         : ClinVarVariantID,
                     'clinvar_condition'          : ClinVarCondition,
                     'clinvar_review_status'      : ClinVarReviewStatus,
                     'clinvar_classification'     : ClinVarClassification,

                    }
