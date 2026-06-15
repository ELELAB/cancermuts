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
from collections import defaultdict


class SequenceProperty(object):
    description = None
    header = None
    category = None
    code = None
    name = None

    def __init__(self):
        self.data_by_pos = defaultdict(list)

    def add_entries(self, positions, sources=None, value=None, **metadata):

        if positions is None:
            positions = []

        if not isinstance(positions, list):
            raise TypeError("positions must be a list of 1-based residue numbers")

        if sources is None:
            sources = []

        entry = {"positions": positions,
                 "sources": sources,
                 "value": value,
                 "metadata": metadata}

        for position in positions:
            self.data_by_pos[position].append(entry)

        return entry

    def get_entries_at(self, position):
        return self.data_by_pos.get(position, [])

    def get_value_at(self,position):
        entries = self.get_entries_at(position)
        if len(entries) == 0:
            return None

        values = []
        for entry in entries:
            value = self.get_value_str(entry)
            values.append(value)
        return values

    def get_sources_at(self, position):
        sources = []
        for entry in self.get_entries_at(position):
            for source in entry.get("sources", []):
                if source not in sources:
                    sources.append(source)
        return sources

    def get_value_str(self, entry):
        return entry.get("value")

    def __getitem__(self, position):
        return self.get_entries_at(position)

    def __repr__(self):
        return "Sequence property '%s' of category '%s' with data at positions %s" % (self.description, self.category,
                                                                                      sorted(self.data_by_pos.keys()))
    def __str__(self):
        return "Sequence property '%s'" % self.description


class LinearMotif(SequenceProperty):
    description = "Residue part of a linear motif"
    header = "linear_motif"
    category = 'linear_motif'

    def add_entries(self, positions, sources=None, name="", id="", **metadata):
        metadata["type"] = name
        metadata["id"] = id
        return super(LinearMotif, self).add_entries(positions=positions, sources=sources, **metadata)

    def get_value_str(self, entry):
        return "%s (%s), %s-%s, %s" % (entry["metadata"]["type"], entry["metadata"]["id"], entry["positions"][0],
                                       entry["positions"][-1], ",".join(s.name for s in entry["sources"]))

class Structure(SequenceProperty):
    description = "Structure"
    header = "structure"
    category = 'structure'

    def add_entries(self, positions, sources=None, name="", **metadata):
        metadata["type"] = name
        return super(Structure, self).add_entries(positions=positions,sources=sources,**metadata)

    def get_value_str(self, entry):
        return "%s, %d-%d, %s" % (entry["metadata"]["type"], entry["positions"][0],
                                  entry["positions"][-1], ",".join(s.name for s in entry["sources"]))

class PTMSite(SequenceProperty):
    def get_value_str(self, entry):
        return self.code

class PhosphorylationSite(PTMSite):
    description = "Phosphorylation site"
    header = "phosphorylation_site"
    category = "ptm_phosphorylation"
    code = "P"

class MethylationSite(PTMSite):
    description = "Methylation site"
    header = "methylation_site"
    category = "ptm_methylation"
    code = "Me"

class AcetylationSite(PTMSite):
    description = "Acetylation site"
    header = "acetylation_site"
    category='ptm_acetylation'
    code = 'Ac'

class SNitrosylationSite(PTMSite):
    description = "S-Nitrosylation site"
    header = "s-nitrosylation_site"
    category='ptm_nitrosylation'
    code = "SN"

class GlycosylationSite(PTMSite):
    description = "Glycosylation site"
    header = "glycosylation_site"
    category='ptm_glycosylation'
    code = "Gly"

    def add_entries(self, positions, sources=None, subtype=None, **metadata):
        metadata["subtypes"] = []

        if subtype is not None:
            metadata["subtypes"].append(subtype)

        return super(GlycosylationSite, self).add_entries(positions=positions, sources=sources, **metadata)

class SumoylationSite(PTMSite):
    description = "Sumoylation site"
    header = "sumoyylation_site"
    category='ptm_sumoylation'
    code = "Sumo"

class UbiquitinationSite(PTMSite):
    description = "Ubiquitination site"
    header = "ubiquitination_site"
    category='ptm_ubiquitination'
    code = 'Ubq'

class CleavageSite(PTMSite):
    description = "Caspase cleavage site"
    header = "cleavage_site"
    category='ptm_cleavage'
    code = 'C'

class DisorderPropensity(SequenceProperty):
    description = "Structural disorder"
    header = "disorder_propensity"
    category = "mobidb_disorder_propensity"

    def add_entries(self, positions, sources=None, disorder_state=None, **metadata):
        return super(DisorderPropensity, self).add_entries(positions=positions, value=disorder_state,
                                                   sources=sources, **metadata)
    def get_value_str(self, entry):
        return entry.get("value")

class Structured(DisorderPropensity):
    def add_entries(self, positions, sources=None, **metadata):
        return super(Structured, self).add_entries(positions=positions, sources=sources, disorder_state="S", **metadata)

class Disordered(DisorderPropensity):
    def add_entries(self, positions, sources=None, **metadata):
        return super(Disordered, self).add_entries(positions=positions, sources=sources, disorder_state="D", **metadata)

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
