# core.py - core classes for the cancermuts package
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
Core cancermuts classes --- :mod:`cancermuts.core`
================================================================
core classes that define the basic framework to handle a sequence
and its mutations.

"""

import csv
import logging
from .metadata import *
from .properties import *
from .log import logger_init

class Sequence(object):

    @logger_init
    def __init__(self, gene_id, sequence, source, aliases=None):
        self.gene_id = gene_id
        if aliases is None:
            self.aliases = {}
        else:
            self.aliases = aliases
        self.source = source
        self.sequence = sequence
        self.positions = []
        self.sequence_numbering = []
        self.properties = {}
        for i,p in enumerate(self.sequence):
            self.positions.append(SequencePosition(p, i+1))
            self.sequence_numbering.append(i+1)

    def index2seq(self, idx):
        return self.sequence_numbering[idx]

    def seq2index(self, seqn):
        return self.sequence_numbering.index(seqn)

    def __iter__(self):
        return iter(self.positions)

    def __repr__(self):
        return "<Sequence of %s from %s, %d positions>" % (self.gene_id, self.source.name, len(self.positions))

    def add_property(self, prop):
        if self.properties.has_key(prop.category):
            self.properties[prop.category].append(prop)
            add_type = "appending"
            self.log.debug("adding property %s to sequence of %s (appeding)" % (str(prop), self.gene_id))
        else:
            self.properties[prop.category] = [prop]
            add_type = "new category"
            self.log.debug("adding property %s to sequence of %s (%s)" % (str(prop), self.gene_id, add_type))


class SequencePosition(object):
    @logger_init
    def __init__(self, wt_residue_type, sequence_position, mutations=None, properties=None):
        self.wt_residue_type = wt_residue_type
        self.sequence_position = sequence_position
        if mutations is None:
            self.mutations = []
        if properties is None:
            self.properties = {} # static properties for position (dependent of WT residue)
    
    def add_mutation(self, mut):
        if mut not in self.mutations:

            #print "appendo da capo"
            self.mutations.append(mut)
            #print self.mutations[-1].metadata
            self.log.info("Adding mutation %s to position %s" % (str(mut), self.__repr__()))
        else:
            self.log.info("Mutation %s already in position %s; will just add sources and metadata" % (str(mut), self.__repr__()))
            pos = self.mutations.index(mut)
            self.mutations[pos].sources.extend(mut.sources)
            for k in mut.metadata.keys():
                if k in self.mutations[pos].metadata.keys():
                    self.mutations[pos].metadata[k].extend(mut.metadata[k])
                    self.log.debug("    metadata %s was extended" % k)
                else:
                    self.mutations[pos].metadata[k] = mut.metadata[k]
                    self.log.debug("    metadata %s was added anew" % k)

    def add_property(self, prop):
        if prop.category in self.properties.keys():
            self.log.info("property %s was replaced with %s" % (self.properties[prop.category], prop))
        else:
            self.log.info("added property %s" % str(prop))

        self.properties[prop.category] = prop

    def __repr__(self):
        return "<SequencePosition, residue %s at position %d>" % (self.wt_residue_type, self.sequence_position)


class Mutation(object):
    @logger_init
    def __init__(self, sequence_position, mutated_residue_type, sources=None, metadata=None):

        self.sequence_position = sequence_position
        self.mutated_residue_type = mutated_residue_type
        if sources is None:
            self.sources = []
        else:
            self.sources = sources
        if metadata is None:
            self.metadata = {}
        else:
            self.metadata = metadata

    def add_source(self, source):
        self.sources.append(source)
    def __eq__(self, other):
        return self.sequence_position == other.sequence_position and self.mutated_residue_type == other.mutated_residue_type
    def __repr__(self):
        return "<Mutation %s%d%s from %s>" % (self.sequence_position.wt_residue_type, 
                                            self.sequence_position.sequence_position, 
                                            self.mutated_residue_type,
                                            ",".join([s.name for s in self.sources]))
    def __str__(self):
        return "%s%d%s" % ( self.sequence_position.wt_residue_type, 
                            self.sequence_position.sequence_position, 
                            self.mutated_residue_type )


class Source(object):
    def __init__(self):
        self.name = name
        self.version = version
    def get_sequence(self, gene_id):
        return None
    def get_mutations(self, ene_id):
        return None
    def get_position_properties(self, position):
        return None
    def get_mutation_properties(self, mutation):
        return None

    def __hash__(self):
        return hash((self.name, self.version))

class MetaTable(object):
    def __init__(self, sequence):
        self.sequence = sequence

    def write_csv(self, outfile="metatable.csv", 
                        mutation_metadata=["cancer_study", "cancer_type", "genomic_coordinates", "genomic_mutations", "revel_score"], 
                        position_properties=['phosphorylation','methylation','ubiquitination','cleavage', 's-nitrosylation','acetylation', 'sumoylation'],
                        sequence_properties=['linear_motif']):
        header =  ['Position',]
        header += [position_properties_classes[p].description for p in position_properties]
        sequence_properties_cols_start = len(header)

        header += [sequence_properties_classes[p].description for p in sequence_properties]
        sequence_properties_col = range(sequence_properties_cols_start, len(header))
        header += ['WT residue', 'Mutated residue', 'Sources']
        for md in mutation_metadata:
            header.append(metadata_classes[md].description)

        with open(outfile, 'wb') as csvfile:

            csvw = csv.writer(csvfile, dialect='excel', delimiter=';')
            csvw.writerow(header)
            csv_rows = []
            positions_mutlist = []

            for gi,p in enumerate(self.sequence.positions):
                base_row = [p.sequence_position]
                val = '-'
                for r in position_properties:
                    try:
                        val = p.properties[r].get_value_str()
                    except:
                        val = '-'
                    base_row.append(val)
                base_row.extend(['-']*len(sequence_properties_col))
                base_row.append(p.wt_residue_type)

                mut_strings = [str(m) for m in p.mutations]
                mut_strings_order = sorted(range(len(mut_strings)), key=mut_strings.__getitem__)
                print mut_strings_order
                for m in mut_strings_order:
                    this_row = list(base_row)
                    this_row.append(p.mutations[m].mutated_residue_type)
                    this_row.append(",".join([s.name for s in p.mutations[m].sources]))
                    for md in mutation_metadata:
                        try:
                            md_values = []
                            if mutation_metadata == 'genomic_mutations':
                                print p.mutations[m].metadata[md], 'UUU'
                            for single_md in p.mutations[m].metadata[md]:
                                md_values.append(single_md.get_value_str())

                            md_str = ", ".join(sorted(list(set(md_values))))
                        except:
                            md_str = "-"
                        this_row.append(md_str)
                    csv_rows.append(this_row)
                    #print "appending", gi
                    positions_mutlist.append(gi)
                    #print positions_mutlist
                if len(mut_strings_order) == 0:
                    this_row = list(base_row)
                    positions_mutlist.append(gi)
                    csv_rows.append(this_row)

            property_rows = {}
            for p in sequence_properties:
                property_rows[p] = [[] for i in self.sequence.positions]

            for pidx, property_name in enumerate(sequence_properties):
                if property_name not in self.sequence.properties.keys():
                    continue    
                for i,p in enumerate(self.sequence.properties[property_name]):
                    if p.category=='linear_motif':
                        sources_str = ",".join(s.name for s in p.sources)
                        cat_str = "%s, %d, %d-%d, %s" % (p.type,
                                                    i,
                                                    p.positions[ 0].sequence_position,
                                                    p.positions[-1].sequence_position,
                                                    sources_str)
                        positions = p.positions

                    for pos in positions:
                        property_rows[property_name][self.sequence.seq2index(pos.sequence_position)].append(cat_str)

                for i,p in enumerate(property_rows[property_name]):
                    indices = [j for j, x in enumerate(positions_mutlist) if x == i]
                    print indices
                    property_str= "|".join(p)
                    for idx in indices:
                        csv_rows[idx][sequence_properties_col[pidx]] = property_str

            for r in csv_rows:
                csvw.writerow(r)
