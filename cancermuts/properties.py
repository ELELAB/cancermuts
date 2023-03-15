    # metadata.py - properties handling for the cancermuts package
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
Classes to handle properties of sequence or positions

"""

class SequenceProperty(object):
    def __init__(self, name, positions=None, sources=None, values=None, metadata=None):
        self.name = name

        if positions is None:
            self.positions = []
        else:
            self.positions = positions

        if sources is None:
            self.sources = []
        else:
            self.sources = sources
        if values is None:
            self.values = []
        else:
            self.values = values
        if metadata is None:
            self.metadata = {}
        else:
            self.metadata = metadata

    def __repr__(self):
        sources_str = ", ".join([x.name for x in self.sources])
        return "<SequenceProperty %s from %s, positions %s>" % (self.name, sources_str, ",".join([str(s.sequence_position) for s in self.positions]))

class PositionProperty(object):
    def __init__(self, name, position, sources=None, values=None, metadata=None):
        self.name = name
        self.position = position
        if sources is None:
            self.sources = []
        else:
            self.sources = sources
        if values is None:
            self.values = []
        else:
            self.values = values
        if metadata is None:
            self.metadata = {}
        else:
            self.metadata = metadata

    def __repr__(self):
        sources_str = ", ".join([x.name for x in self.sources])
        return "<PositionProperty %s from %s>" % (self.name, sources_str)

class LinearMotif(SequenceProperty):
    description = "Residue part of a linear motif"
    header = "linear_motif"
    category = 'linear_motif'

    def __init__(self, positions, sources, name="", id=""):
        super(LinearMotif, self).__init__(  name="Linear motif",
                                                    positions=positions,
                                                    sources=sources,
                                                    values=None,
                                                    metadata=None  )
        self.type = name
        self.id = id

    def get_value_str(self):
        return f'{self.type} ({self.id}), {self.positions[ 0].sequence_position}-{self.positions[-1].sequence_position}, {",".join(s.name for s in self.sources)}'


class Structure(SequenceProperty):
    description = "Structure"
    header = "structure"
    category = 'structure'

    def __init__(self, positions, sources, name=""):
        super(Structure, self).__init__(  name="Structure",
                                                    positions=positions,
                                                    sources=sources,
                                                    values=None,
                                                    metadata=None  )
        self.type = name

    def get_value_str(self):
        return "%s, %d-%d, %s" % (    self.type,
                                      self.positions[ 0].sequence_position,
                                      self.positions[-1].sequence_position,
                                      ",".join(s.name for s in self.sources))



class PhosphorylationSite(PositionProperty):
    description = "Phosphorylation site"
    header = "phosphorylation_site"
    category="ptm_phosphorylation"
    code = "P"

    def __init__(self, position, sources):
        super(PhosphorylationSite, self).__init__(  name="Phosphorylation Site",
                                                    position=position,
                                                    sources=sources,
                                                    values={},
                                                    metadata={}  )

    def get_value_str(self):
        return self.code

class MethylationSite(PositionProperty):
    description = "Methylation site"
    header = "methylation_site"
    category="ptm_methylation"
    code = "Me"

    def __init__(self, position, sources):
        super(MethylationSite, self).__init__(  name="Methylation Site",
                                                    position=position,
                                                    sources=sources,
                                                    values={},
                                                    metadata={}  )

    def get_value_str(self):
        return self.code

class AcetylationSite(PositionProperty):
    description = "Acetylation site"
    header = "acetylation_site"
    category='ptm_acetylation'
    code = 'Ac'

    def __init__(self, position, sources):
        super(AcetylationSite, self).__init__(  name="Acetylation Site",
                                                    position=position,
                                                    sources=sources,
                                                    values={},
                                                    metadata={}  )

    def get_value_str(self):
        return self.code

class SNitrosylationSite(PositionProperty):
    description = "S-Nitrosylation site"
    header = "s-nitrosylation_site"
    category='ptm_nitrosylation'
    code = "SN"

    def __init__(self, position, sources):
        super(SNitrosylationSite, self).__init__(  name="S-Nytrosilation site",
                                                    position=position,
                                                    sources=sources,
                                                    values={},
                                                    metadata={}  )

    def get_value_str(self):
        return self.code

class OGalNAcSite(PositionProperty):
    description = "OGalNAc site"
    header = "ogalnac_site"
    category='ptm_ogalnac'
    code = "O-GalNAc"

    def __init__(self, position, sources):
        super(OGalNAcSite, self).__init__(  name="o-GalNAc Site",
                                                    position=position,
                                                    sources=sources,
                                                    values={},
                                                    metadata={}  )

    def get_value_str(self):
        return self.code

class OGlcNAcSite(PositionProperty):
    description = "OGlcNAc site"
    header = "oglcnac_site"
    category='ptm_oglcnac'
    code = "O-GlcNAc"


    def __init__(self, position, sources):
        super(OGlcNAcSite, self).__init__(  name="o-GlcNAc Site",
                                                    position=position,
                                                    sources=sources,
                                                    values={},
                                                    metadata={}  )

    def get_value_str(self):
        return self.code

class SumoylationSite(PositionProperty):
    description = "Sumoylation site"
    header = "sumoyylation_site"
    category='ptm_sumoylation'
    code = "Sumo"


    def __init__(self, position, sources):
        super(SumoylationSite, self).__init__(  name="Sumoylation Site",
                                                    position=position,
                                                    sources=sources,
                                                    values={},
                                                    metadata={}  )
    def get_value_str(self):
        return self.code

class UbiquitinationSite(PositionProperty):
    description = "Ubiquitination site"
    header = "ubiquitination_site"
    category='ptm_ubiquitination'
    code = 'Ubq'

    def __init__(self, position, sources):
        super(UbiquitinationSite, self).__init__(  name="Ubiquitination Site",
                                                    position=position,
                                                    sources=sources,
                                                    values={},
                                                    metadata={}  )
    def get_value_str(self):
        return self.code

class CleavageSite(PositionProperty):
    description = "Caspase cleavage site"
    header = "cleavage_site"
    category='ptm_cleavage'
    code = 'C'

    def __init__(self, position, sources):
        super(CleavageSite, self).__init__(  name="Caspase cleavage site",
                                                    position=position,
                                                    sources=sources,
                                                    values={},
                                                    metadata={}  )

    def get_value_str(self):
        return self.code


class DisorderPropensity(PositionProperty):
    description = "Structural disorder"
    header = "disorder_propensity"
    category = 'disorder_propensity'

#   def __init__(self, name, category, position, sources=None, values=None, metadata=None):
    def __init__(self, position, sources, disorder_state):
        super(DisorderPropensity, self).__init__(name="Disorder Propensity",
                                                 position=position,
                                                 sources=sources,
                                                 values={},
                                                 metadata={})
        self.disorder_propensity = disorder_state

    def get_value(self):
        return self.disorder_propensity

    def get_value_str(self):
        return self.disorder_propensity

    def __repr__(self):
        return "<StructuralDisorder, %s>" % self.disorder_propensity

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        return self.source == other.source and \
               self.disorder_propensity == other.disorder_propensity

    def __hash__(self):
        return hash((self.source, self.disorder_propensity))

class Structured(DisorderPropensity):
    def __init__(self, position, sources):
        super(Ordered, self).__init__(name="Disorder Propensity",
                                                 position=position,
                                                 sources=sources,
                                                 disorder_state='S')

class Disordered(DisorderPropensity):
    def __init__(self, position, sources):
        super(Ordered, self).__init__(name="Disorder Propensity",
                                                 position=position,
                                                 sources=sources,
                                                 disorder_state='D')

position_properties_classes = {  'ptm_phosphorylation'            : PhosphorylationSite,
                                 'ptm_methylation'                : MethylationSite,
                                 'ptm_acetylation'                : AcetylationSite,
				                 'ptm_nitrosylation'              : SNitrosylationSite,
                                 'ptm_ogalnac'                    : OGalNAcSite,
                                 'ptm_oglcnac'                    : OGlcNAcSite,
                                 'ptm_sumoylation'                : SumoylationSite,
                                 'ptm_ubiquitination'             : UbiquitinationSite,
                                 'ptm_cleavage'                   : CleavageSite,
                                 'mobidb_disorder_propensity'     : DisorderPropensity
                              }

sequence_properties_classes = {  'linear_motif'                   : LinearMotif,
                                 'structure'            : Structure
                              }