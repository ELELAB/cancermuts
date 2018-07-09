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

    def __init__(self, source, genome_version, chromosome, coord, wt, mut):
        super(GenomicMutation, self).__init__(source)
        if mut not in self.allowed_bases or wt not in self.allowed_bases:
            raise TypeError
        self.genome_version = genome_version
        self.chr = chromosome
        self.coord = coord
        self.wt = wt
        self.mut = mut

    def get_value(self):
        return [self.genome_version, self.chr, self.coord, self.wt, self.mut]

    def get_value_str(self):
        return "%s,chr%s:%s%s>%s" % (self.genome_version, self.chr, self.coord, self.wt, self.mut)

    def __repr__(self):
        return "<GenomicMutation %s from %s>" % (self.get_value_str(), self.source.name)

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        return self.source == other.source and \
        self.chr == other.chr and \
        self.coord == other.coord and \
        self.wt == other.wt and \
        self.mut == other.mut

    def __hash__(self):
        return hash((self.source, self.chr, self.coord, self.wt, self.mut))

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

metadata_classes = {'cancer_type'         : CancerType,
                    'cancer_study'        : CancerStudy,
                    'genomic_coordinates' : GenomicCoordinates,
                    'genomic_mutations'   : GenomicMutation,
                    'revel_score'         : DbnsfpRevel}