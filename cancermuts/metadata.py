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

    def __init__(self, source, genome_version, chromosome, coord_start, coord_end, ref):
        super(GenomicCoordinates, self).__init__(source)
        self.genome_version = genome_version
        self.chr = chromosome
        self.coord_start = coord_start
        self.coord_end = coord_end
        self.ref = ref

    def get_value(self):
        return [self.genome_version, self.chr, self.coord_start, self.coord_end, self.ref]

    def get_value_str(self):
        return "%s,chr%s:%s-%s" % (self.genome_version, self.chr, self.coord_start, self.coord_end)

    def get_coord(self):
        return self.coord_start

    def __repr__(self):
        return "<GenomicCoordinates %s:%s-%s in %s from %s>" % (self.chr, self.coord_start, self.coord_end, self.genome_version, self.source.name)

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        return self.source == other.source and \
        self.chr == other.chr and \
        self.coord_start == other.coord_start and \
        self.coord_end == other.coord_end and \
        self.ref == other.ref and \
        self.genome_version == other.genome_version

    def __hash__(self):
        return hash((self.source, self.chr, self.coord_start, self.coord_end, self.ref, self.genome_version))

class GenomicMutation(Metadata):

    description = "Genomic mutation"

    allowed_bases = ['A', 'C', 'G', 'T']

    complementarity = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

    @logger_init
    def __init__(self, source, genome_version, chromosome, strand, coord, wt, mut):
        super(GenomicMutation, self).__init__(source)
        if mut not in self.allowed_bases or wt not in self.allowed_bases:
            raise TypeError
        self.genome_version = genome_version
        self.chr = chromosome
        self.coord = coord
        self.wt = wt
        self.mut = mut
        if strand == '':
            self.log.info("strand will be defaulted to + for %s" % self)
            strand = '+'
        self.strand = strand

    def get_value(self):
        return [self.genome_version, self.chr, self.coord, self.strand, self.wt, self.mut]

    def get_value_str(self, fmt='csv'):
        if fmt == 'csv':
            return "%s,chr%s:%s%s>%s" % (self.genome_version, self.chr, self.coord, self.get_sense_wt(), self.get_sense_mut())
        if fmt == 'gnomad':
            return '%s-%s-%s-%s' % (self.chr, self.coord, self.get_sense_wt(), self.get_sense_mut())
        return None

    def get_coord(self):
        return self.coord

    def get_sense_wt(self):
        return self._get_sense_base(self.wt)

    def get_sense_mut(self):
        return self._get_sense_base(self.mut)

    def _get_sense_base(self, base):
        if self.strand == '+':
            return base
        elif self.strand == '-':
            return self.complementarity[base]

    def as_hg19(self):
        if self.genome_version == 'hg19':
            return GenomicMutation(self.source, self.genome_version, self.chr, self.strand, self.coord, self.wt, self.mut)
        elif self.genome_version == 'hg38':
            converted_coords = lo_hg38_hg19.convert_coordinate('chr%s' % self.chr, int(self.get_coord()))
            assert len(converted_coords) == 1
            return GenomicMutation(self.source, 'hg19', converted_coords[0][0][3:], self.strand, converted_coords[0][1], self.wt, self.mut)
        else:
            raise TypeError

    def as_hg38(self):
        if self.genome_version == 'hg38':
            return GenomicMutation(self.source, self.genome_version, self.chr, self.strand, self.coord, self.wt, self.mut)
        elif self.genome_version == 'hg19':
            converted_coords = lo_hg19_hg38.convert_coordinate('chr%s' % self.chr, int(self.get_coord()))
            assert len(converted_coords) == 1
            return GenomicMutation(self.source, 'hg38', converted_coords[0][0][3:], self.strand, converted_coords[0][1], self.wt, self.mut)
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
        self.chr == other.chr and \
        self.strand == other.strand and \
        self.coord == other.coord and \
        self.wt == other.wt and \
        self.mut == other.mut

    def __hash__(self):
        return hash((self.source, self.chr, self.strand, self.coord, self.wt, self.mut))

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

    @classmethod
    def set_version_in_desc(cls, version_string):
        cls.description = "%s (%s)" % (cls.description, version_string)

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
    description = "%s allele frequency" % freq_type

    def __init__(self, source, frequency):
        super(gnomADExomeAlleleFrequency, self).__init__(source, frequency)

    def __repr__(self):
        return "<gnomADExomeAlleleFrequency, %f>" % self.frequency

class gnomADGenomeAlleleFrequency(gnomADAlleleFrequency):

    freq_type = "Genome"
    description = "%s allele frequency" % freq_type

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
