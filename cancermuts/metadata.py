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

    def __init__(self, source, genome_build, chromosome, coord_start, coord_end, ref):
        super(GenomicCoordinates, self).__init__(source)
        self.genome_build = genome_build
        self.chr = chromosome
        self.coord_start = coord_start
        self.coord_end = coord_end
        self.ref = ref

    def get_value(self):
        return [self.genome_build, self.chr, self.coord_start, self.coord_end, self.ref]

    def get_value_str(self):
        return "%s,chr%s:%s-%s" % (self.genome_build, self.chr, self.coord_start, self.coord_end)

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

    allowed_bases = ['A', 'C', 'G', 'T']

    _mut_snv_regexp = '^[0-9]+:g\.[0-9]+[ACTG]>[ACTG]'
    _mut_snv_prog = re.compile(_mut_snv_regexp)
    _mut_snv_parse = '{chr:d}:g.{coord:d}{ref:l}>{alt:l}'

    @logger_init
    def __init__(self, source, genome_build, definition):
        super(GenomicMutation, self).__init__(source)

        self.genome_build = genome_build

        if self._mut_snv_prog.match(definition):
            tokens = parse(self._mut_snv_parse, definition)
            if tokens['ref'] not in self.allowed_bases or \
               tokens['alt'] not in self.allowed_bases:
                self.log.warning(f'this mutation does not specify allowed nucleotides:  {genome_build}, {chromosome}, {strand}, {coord}, {wt}, {mut}')
                return None

            self.chr = tokens['chr']
            self.coord = tokens['coord']
            self.ref = tokens['ref']
            self.alt = tokens['alt']
            self.definition = definition
            self.is_snv = True

        else:
            self.chr = None
            self.coord = None
            self.ref = None
            self.alt = None
            self.definition = definition
            self.is_snv = False

    def get_value(self):
        return [self.genome_build, self.chr, self.coord, self.wt, self.mut]

    def get_value_str(self, fmt='csv'):
        if fmt == 'csv':
            return f"{self.genome_build},{self.definition}"
        if fmt == 'gnomad':
            if self.is_snv:
                return f"{self.chr}-{self.coord}-{self.wt}-{self.mut}"
            else:
                return None
        else:
            return None

    def get_coord(self):
        return self.coord

    def as_hg19(self):
        if self.is_snv:
            if self.genome_build == 'hg19':
                return GenomicMutation(self.source, self.genome_build, self.chr, self.coord, self.wt, self.mut)
            elif self.genome_build == 'hg38':
                converted_coords = lo_hg38_hg19.convert_coordinate('chr%s' % self.chr, int(self.get_coord()))
                assert len(converted_coords) == 1
                return GenomicMutation(self.source, 'hg19', converted_coords[0][0][3:], converted_coords[0][1], self.wt, self.mut)
            else:
                raise TypeError

    def as_hg38(self):
        if self.is_snv:
            if self.genome_build == 'hg38':
                return GenomicMutation(self.source, self.genome_build, self.chr, self.coord, self.wt, self.mut)
            elif self.genome_build == 'hg19':
                converted_coords = lo_hg19_hg38.convert_coordinate('chr%s' % self.chr, int(self.get_coord()))
                assert len(converted_coords) == 1
                return GenomicMutation(self.source, 'hg38', converted_coords[0][0][3:], converted_coords[0][1], self.wt, self.mut)
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

class DbnsfpRevel(Metadata):

    description = "Revel score"

    def __init__(self, source, score):
        super(DbnsfpRevel, self).__init__(source)
        self.source = source
        self.score = score

    def get_value(self):
        return self.score

    def get_value_str(self):
        return "%.3f" % self.score

    def __repr__(self):
        return "<DbnsfpRevel, %.3f>" % self.score

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

    @classmethod
    def set_version_in_desc(cls, version_string):
        cls.description = "%s (%s)" % (cls.basic_description, version_string)

    def __init__(self, source, frequency):
        super(gnomADGenomeAlleleFrequency, self).__init__(source, frequency)

    def __repr__(self):
        return "<gnomADGenomeAlleleFrequency, %f>" % self.frequency

class CancerSite(Metadata):

    description = "Cancer site"
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

metadata_classes = { 
                     'cancer_type'                 : CancerType,
                     'cancer_study'                : CancerStudy,
                     'genomic_coordinates'         : GenomicCoordinates,
                     'genomic_mutations'           : GenomicMutation,
                     'revel_score'                 : DbnsfpRevel,
                     'gnomad_genome_allele_frequency' : gnomADGenomeAlleleFrequency,
                     'gnomad_exome_allele_frequency' : gnomADExomeAlleleFrequency,
                     'cancer_site'                 : CancerSite,
                     'cancer_histology'            : CancerHistology,
                   }
