    # metadata.py - properties handling for the cancermuts package
# (c) 2018 Matteo Tiberti <matteo.tiberti@gmail.com>
# (c) 2023 Katrine Meldgård <katrine@meldgaard.dk>
# (c) 2026 Beatrice Drago
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

        if not isinstance(self.positions, list):
            raise TypeError("positions must be a list of 1-based residue numbers")

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
        positions_str = ",".join([str(pos) for pos in self.positions])

        return "<SequenceProperty %s from %s, positions %s>" % (self.name, sources_str, positions_str)


class LinearMotif(SequenceProperty):
    description = "Residue part of a linear motif"
    header = "linear_motif"
    category = 'linear_motif'

    def __init__(self, positions, sources, name="", id=""):
        super(LinearMotif, self).__init__(  name="Linear motif",
                                                    positions=positions,
                                                    sources=sources,
                                                    values=None,
                                                    metadata=None )
        self.type = name
        self.id = id

    def get_value_str(self):
        return f'{self.type} ({self.id}), {self.positions[ 0]}-{self.positions[-1]}, {",".join(s.name for s in self.sources)}'


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
                                      self.positions[ 0],
                                      self.positions[-1],
                                      ",".join(s.name for s in self.sources))



class PhosphorylationSite(SequenceProperty):
    description = "Phosphorylation site"
    header = "phosphorylation_site"
    category="ptm_phosphorylation"
    code = "P"

    def __init__(self, positions, sources):
        super(PhosphorylationSite, self).__init__(  name="Phosphorylation Site",
                                                    positions=positions,
                                                    sources=sources,
                                                    values={},
                                                    metadata={}  )

    def get_value_str(self):
        return self.code

class MethylationSite(SequenceProperty):
    description = "Methylation site"
    header = "methylation_site"
    category="ptm_methylation"
    code = "Me"

    def __init__(self, positions, sources):
        super(MethylationSite, self).__init__(  name="Methylation Site",
                                                    positions=positions,
                                                    sources=sources,
                                                    values={},
                                                    metadata={}  )

    def get_value_str(self):
        return self.code

class AcetylationSite(SequenceProperty):
    description = "Acetylation site"
    header = "acetylation_site"
    category='ptm_acetylation'
    code = 'Ac'

    def __init__(self, positions, sources):
        super(AcetylationSite, self).__init__(  name="Acetylation Site",
                                                    positions=positions,
                                                    sources=sources,
                                                    values={},
                                                    metadata={}  )

    def get_value_str(self):
        return self.code

class SNitrosylationSite(SequenceProperty):
    description = "S-Nitrosylation site"
    header = "s-nitrosylation_site"
    category='ptm_nitrosylation'
    code = "SN"

    def __init__(self, positions, sources):
        super(SNitrosylationSite, self).__init__(  name="S-Nytrosilation site",
                                                    positions=positions,
                                                    sources=sources,
                                                    values={},
                                                    metadata={}  )

    def get_value_str(self):
        return self.code

class GlycosylationSite(SequenceProperty):
    description = "Glycosylation site"
    header = "glycosylation_site"
    category='ptm_glycosylation'
    code = "Gly"

    def __init__(self, positions, sources):
        super(GlycosylationSite, self).__init__(  name="Glycosylation Site",
                                                    positions=positions,
                                                    sources=sources,
                                                    values={},
                                                    metadata={"subtypes": []}  )

    def add_subtype(self,subtype):
        if subtype not in self.metadata["subtypes"]:
            self.metadata["subtypes"].append(subtype)

    def get_value_str(self):
        return self.code

class SumoylationSite(SequenceProperty):
    description = "Sumoylation site"
    header = "sumoyylation_site"
    category='ptm_sumoylation'
    code = "Sumo"

    def __init__(self, positions, sources):
        super(SumoylationSite, self).__init__(  name="Sumoylation Site",
                                                    positions=positions,
                                                    sources=sources,
                                                    values={},
                                                    metadata={}  )
    def get_value_str(self):
        return self.code

class UbiquitinationSite(SequenceProperty):
    description = "Ubiquitination site"
    header = "ubiquitination_site"
    category='ptm_ubiquitination'
    code = 'Ubq'

    def __init__(self, positions, sources):
        super(UbiquitinationSite, self).__init__(  name="Ubiquitination Site",
                                                    positions=positions,
                                                    sources=sources,
                                                    values={},
                                                    metadata={}  )
    def get_value_str(self):
        return self.code

class CleavageSite(SequenceProperty):
    description = "Caspase cleavage site"
    header = "cleavage_site"
    category='ptm_cleavage'
    code = 'C'

    def __init__(self, positions, sources):
        super(CleavageSite, self).__init__(  name="Caspase cleavage site",
                                                    positions=positions,
                                                    sources=sources,
                                                    values={},
                                                    metadata={}  )

    def get_value_str(self):
        return self.code

class DisorderPropensity(SequenceProperty):
    description = "Structural disorder"
    header = "disorder_propensity"
    category = 'disorder_propensity'

    def __init__(self, positions, sources, disorder_state):
        super(DisorderPropensity, self).__init__(name="Disorder Propensity",
                                                 positions=positions,
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
        return self.sources == other.sources and \
               self.disorder_propensity == other.disorder_propensity

    def __hash__(self):
        return hash((self.sources, self.disorder_propensity))

class Structured(DisorderPropensity):
    def __init__(self, positions, sources):
        super(Ordered, self).__init__(name="Disorder Propensity",
                                                 positions=positions,
                                                 sources=sources,
                                                 disorder_state='S')

class Disordered(DisorderPropensity):
    def __init__(self, positions, sources):
        super(Ordered, self).__init__(name="Disorder Propensity",
                                                 positions=positions,
                                                 sources=sources,
                                                 disorder_state='D')

sequence_properties_classes = {  'ptm_phosphorylation'            : PhosphorylationSite,
                                 'ptm_methylation'                : MethylationSite,
                                 'ptm_acetylation'                : AcetylationSite,
				                 'ptm_nitrosylation'              : SNitrosylationSite,
                                 'ptm_glycosylation'              : GlycosylationSite,
                                 'ptm_sumoylation'                : SumoylationSite,
                                 'ptm_ubiquitination'             : UbiquitinationSite,
                                 'ptm_cleavage'                   : CleavageSite,
                                 'mobidb_disorder_propensity'     : DisorderPropensity,
                                 'linear_motif'                   : LinearMotif,
                                 'structure'                      : Structure
                              }
