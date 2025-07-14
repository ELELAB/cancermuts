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
    """The most fundamental class of Cancermuts, this class starts from a
    protein sequence definition. It acts as a collection of ordered 
    `SequencePosition` objects which depend on the specific
    sequence itself.

    Parameters
    ----------
    gene_id : :obj:`str`
        ID of the gene (gene name) to which the sequence belongs
    sequence : `str`
        Protein sequence corresponding to the gene and isoform of interest
    source: :obj:`cancermuts.datasources.DataSource`
        source from where the Sequence has been downloaded from
    aliases : :obj:`dict` of :obj:`str`, both for key and value, or None, optional
        This argument assigns the `aliases` dictionary for the Sequence class.
        if None, an empty dictionary is created. 

    Attributes
    ----------
    source : :obj:`cancermuts.datasources.Datasource`
        Source from where the protein sequence is downloaded from
    gene_id : :obj:`int`, optional
        gene name to which the sequence to be downloaded belongs to
    sequence : :obj:`str`
        protein sequence as obtained from the data source
    sequence_numbering : :obj:`list of int`
        list of integers, starting from one, each representing a position
    properties : :obj:`dict`
        dictionary including the downloaded protein-associated properties.
        These can span one or more residues.
    aliases : :obj:`dict`
        This dictionary maps different "aliases" (i.e. identifiers) for the
        protein of interest. The key should represent the type of identifier
        it's being stored in the value and can be any string. For instance,
        one could have "cosmic" : "AMBRA1" to mean that the protein AMBRA1
        has identifier AMBRA1 in cosmic
    """

    @logger_init
    def __init__(self, gene_id, uniprot_ac, sequence, source, isoform=None, is_canonical=False, aliases=None):
        self.gene_id = gene_id
        self.uniprot_ac = uniprot_ac
        self.isoform = isoform
        self.is_canonical = is_canonical

        if aliases is None:
            self.aliases = {}
        else:
            if "uniprot" in aliases and aliases["uniprot"] != uniprot_ac:
                raise TypeError(
                    f"Mismatch between provided uniprot_ac ('{uniprot_ac}') "
                    f"and aliases['uniprot'] ('{aliases['uniprot']}')"
                )
            self.aliases = aliases
            
        self.aliases["uniprot"] = self.uniprot_ac
        self.source = source
        self.sequence = sequence
        self.positions = []
        self.sequence_numbering = []
        self.properties = {}
        for i,p in enumerate(self.sequence):
            self.positions.append(SequencePosition(p, i+1))
            self.sequence_numbering.append(i+1)

    def index2seq(self, idx):
        """
        Returns the sequence numbering corresponding to the `idx`-th
        residue, starting from 0. usually this corresponds to idx+1.

        Parameters
        ----------
        idx : :obj:`int`
            0-index position in the protein sequence

        Returns
        ----------
        i:obj:`int`
            sequence number corresponding to the idx-th position
        """


    def seq2index(self, seqn):
        """
        Returns the 0-indexed position number in the sequence corresponding
        to sequence number `seqn`. It is usually equal to `seqn-1`..

        Parameters
        ----------
        seqn : :obj:`int`
            residue position number in the protein sequence

        Returns
        ----------
        :obj:`int`
            sequential number corresponding to the position in the sequence
            (starting from 0)
        """
        return self.sequence_numbering.index(seqn)

    def __iter__(self):
        return iter(self.positions)

    def __repr__(self):
        return f"<Sequence gene_id={self.gene_id}, uniprot_ac={self.uniprot_ac}, isoform={self.isoform}, is_canonical={self.is_canonical}, source={self.source.name}, {len(self.positions)} positions>"
    
    def add_property(self, prop):
        """
        Adds sequence property to sequence object. If a property of the same
        category is already present, the property will be added to the same
        category; otherwise the category will be created anew

        Parameters
        ----------
        prop : :obj:`cancermuts.properties.SequenceProperty`
            new property to be added
        """

        if prop.category in self.properties:
            self.properties[prop.category].append(prop)
            add_type = "appending"
            self.log.debug("adding property %s to sequence of %s (appeding)" % (str(prop), self.gene_id))
        else:
            self.properties[prop.category] = [prop]
            add_type = "new category"
            self.log.debug("adding property %s to sequence of %s (%s)" % (str(prop), self.gene_id, add_type))


class SequencePosition(object):
    """This class describe a certain sequence position in a protein.
    SequencePositions usually belong to a Sequence object.
    It is possible to annotate a Sequence position with either a Mutation
    or a PositionProperty.

    Attributes
    ----------
    source : :obj:`cancermuts.datasources.Datasource`
    wt_residue_type : :obj:`str`
        single-letter code wild-type residue for this position
    sequence_position : :obj:`int`
        number corresponding to the sequence position for this position
    mutations : :obj:`list` of `cancermuts.core.Mutation` objects
        list of mutations for this sequence position (if any)
    properties : :obj:`dict`
        dictionary encoding position properties. This dictionary needs to have
        :obj:`str` as key and `cancermuts.properties.PositionProperty` as value.
    """
    description = 'Position'
    header = "aa_position"
    @logger_init
    def __init__(self, wt_residue_type, sequence_position, mutations=None, properties=None):
        """Constructor for the SequencePosition class.

        Parameters
        ----------
        wt_residue_type : :obj:`str`
            single-letter code wild-type residue for this position
        sequence_position : :obj:`int`
            number corresponding to the sequence position for this position
        mutations : :obj:`list` of `cancermuts.core.Mutation` objects or None
            list of mutations for this sequence position to be added
        properties : :obj:`dict` or :obj:`None`
            dictionary encoding position properties. This dictionary needs to have
            :obj:`str` as key and `cancermuts.properties.PositionProperty` as value.
        """


        self.wt_residue_type = wt_residue_type
        self.sequence_position = sequence_position
        if mutations is None:
            self.mutations = []
        else:
            self.mutations = mutations
        if properties is None:
            self.properties = {} # static properties for position (dependent of WT residue)
        else:
            self.properties = properties

    def add_mutation(self, mut):
        """
        Adds mutation to a :obj:`SequencePosition` object. If the mutation (in
        terms of amino-acid substitution) is already present, source and 
        metadata will be added to the already-present mutation. Otherwise
        the mutation is added anew.

        Parameters
        ----------
        mut : :obj:`cancermuts.core.Mutation`
            new Mutation object to be added to the SequencePosition
        """

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
        """
        Adds position property to a :obj:`SequencePosition` object. If a
        property of the same category is already present, it will be overwritten.
        Otherwise, the property is just added to the object.

        Parameters
        ----------
        prop : :obj:`cancermuts.properties.PositionProperty`
            new Mutation object to be added to the SequencePosition
        """

        if prop.category in self.properties:
            self.log.info("property %s was replaced with %s" % (self.properties[prop.category], prop))
        else:
            self.log.info("added property %s" % str(prop))

        self.properties[prop.category] = prop

    def __repr__(self):
        return "<SequencePosition, residue %s at position %d>" % (self.wt_residue_type, self.sequence_position)


class Mutation(object):
    """This class describe a missense mutation as amino-acid replacement.

    Attributes
    ----------
    sequence_position : :obj:`cancermuts.core.SequencePosition`
        `SequencePosition` object to which this mutation belongs, corresponding
        to the residue that is mutated
    sources : :obj:`list` of :obj:`cancermuts.datasources.Datasource`
        source the mutation was derived from
    mutated_residue_type : :obj:`str`
        single-letter code for the mutated residue type for this position
    mutations : :obj:`list` of `cancermuts.core.Mutation` objects
        list of mutations for this sequence position (if any)
    properties : :obj:`dict`
        dictionary encoding position properties. This dictionary needs to have
        :obj:`str` as key and `cancermuts.properties.PositionProperty` as value.
    """

    @logger_init
    def __init__(self, sequence_position, mutated_residue_type, sources=None, metadata=None):
        """Constructor for the Mutation class.

        Parameters
        ----------
        sequence_position : :obj:`cancermuts.core.SequencePosition`
            `SequencePosition` object to which this mutation belongs, corresponding
            to the residue that is mutated
            single-letter code wild-type residue for this position
        mutated_residue_type : :obj:`str`
            single-letter code for the mutated residue type for this position
        sources : :obj:`list` of :obj:`cancermuts.datasources.Datasource`
            list of one or more sources the mutation was derived from
        metadata : :obj:`dict` or :obj:`None`
            dictionary encoding position properties. This dictionary needs to have
            :obj:`str` as key and `cancermuts.properties.PositionProperty` as value.
        """

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
