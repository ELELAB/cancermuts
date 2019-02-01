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
        if prop.category in self.properties:
            self.properties[prop.category].append(prop)
            add_type = "appending"
            self.log.debug("adding property %s to sequence of %s (appeding)" % (str(prop), self.gene_id))
        else:
            self.properties[prop.category] = [prop]
            add_type = "new category"
            self.log.debug("adding property %s to sequence of %s (%s)" % (str(prop), self.gene_id, add_type))


class SequencePosition(object):

    description = 'Position'
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
            for k in mut.metadata:
                if k in self.mutations[pos].metadata:
                    self.mutations[pos].metadata[k].extend(mut.metadata[k])
                    self.log.debug("    metadata %s was extended" % k)
                else:
                    self.mutations[pos].metadata[k] = mut.metadata[k]
                    self.log.debug("    metadata %s was added anew" % k)

    def add_property(self, prop):
        if prop.category in self.properties:
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

