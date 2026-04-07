# core.py - core classes for the cancermuts package
# (c) 2018 Matteo Tiberti <matteo.tiberti@gmail.com>
# (c) 2025 Pablo Sanchez-Izquierdo
# This file is part of cancermuts
#
# cancermuts is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cancermuts is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with cancermuts.  If not, see <http://www.gnu.org/licenses/>.


"""
Core cancermuts classes --- :mod:`cancermuts.core`
================================================================
core classes that define the basic framework to handle a sequence
and its variants.

"""

import csv
import logging
from .metadata import *
from .properties import *
from .log import logger_init

class Sequence(object):
    """The most fundamental class of Cancermuts, this class starts from a
    protein sequence definition. It acts as a collection of ordered 
    `SequencePosition` objects and protein variant objects.

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
            if "uniprot_acc" in aliases and aliases["uniprot_acc"] != uniprot_ac:
                raise TypeError(
                    f"Mismatch between provided uniprot_ac ('{uniprot_ac}') "
                    f"and aliases['uniprot_acc'] ('{aliases['uniprot_acc']}')"
                )
            self.aliases = aliases
        
        self.aliases["uniprot_acc"] = self.uniprot_ac
        self.source = source
        self.sequence = sequence
        self.variants = []
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

    def add_variant(self, var):
        """
        Adds variant to a :obj:`Sequence` object. If the variant is already present, 
        source and metadata will be added to the already-present variant. Otherwise
        the variant is added as new.

        Parameters
        ----------
        var : :obj:`cancermuts.core.ProteinVariant`
            new ProteinVariant object to be added to the Sequence
        """

        if var not in self.variants:

            self.variants.append(var)
            self.log.info("Adding variant %s to sequence %s" % (str(var), self.__repr__()))
        else:
            self.log.info("Variant %s already in sequence %s; will just add sources and metadata" % (str(var), self.__repr__()))
            idx = self.variants.index(var)
            self.variants[idx].sources.extend(var.sources)
            for k in var.metadata:
                if k in self.variants[idx].metadata:
                    self.variants[idx].metadata[k].extend(var.metadata[k])
                    self.log.debug("    metadata %s was extended" % k)
                else:
                    self.variants[idx].metadata[k] = var.metadata[k]
                    self.log.debug("    metadata %s was added anew" % k)


class SequencePosition(object):
    """This class describes a certain sequence position in a protein.
    SequencePositions belong to a Sequence object and represent
    the wild-type residue at a given sequence coordinate.

    Attributes
    ----------
    source : :obj:`cancermuts.datasources.Datasource`
    wt_residue_type : :obj:`str`
        single-letter code wild-type residue for this position
    sequence_position : :obj:`int`
        number corresponding to the sequence position for this position
    """
    description = 'Position'
    header = "aa_position"
    @logger_init
    def __init__(self, wt_residue_type, sequence_position):
        """Constructor for the SequencePosition class.

        Parameters
        ----------
        wt_residue_type : :obj:`str`
            single-letter code wild-type residue for this position
        sequence_position : :obj:`int`
            number corresponding to the sequence position for this position
        """

        self.wt_residue_type = wt_residue_type
        self.sequence_position = sequence_position

    def __repr__(self):
        return "<SequencePosition, residue %s at position %d>" % (self.wt_residue_type, self.sequence_position)

class ProteinVariant(object):
    """This class describes in-frame protein variants on a reference sequence.

        Attributes
        ----------
        sequence : :obj:`cancermuts.core.Sequence`
            sequence to which this variant belongs
        start : :obj:`int`
            1-based start coordinate on the protein sequence
        end : :obj:`int`
            1-based end coordinate on the protein sequence
        ref : :obj:`str`
            reference amino-acid sequence
        alt : :obj:`str`
            altered amino-acid sequence
        variant_type: :obj:`str`
            Type of protein variant. Supported values currently include
            ``"substitution"``, ``"deletion"``, ``"insertion"``, and ``"delins"``.
        sources : :obj:`list` of :obj:`cancermuts.datasources.Datasource`
            sources the variant was derived from
        metadata : :obj:`dict`
            Dictionary encoding variant-associated metadata.
    """

    @logger_init
    def __init__(self, sequence, start, end, ref, alt, variant_type, sources=None, metadata=None):
        """Constructor for the ProteinVariant class.
        
        Parameters
        ----------
        sequence : :obj:`cancermuts.core.Sequence`
            sequence to which this variant belongs
        start : :obj:`int`
            1-based start coordinate on the protein sequence
        end : :obj:`int`
            1-based end coordinate on the protein sequence
        ref : :obj:`str`
            reference amino-acid sequence
        alt : :obj:`str`
            altered amino-acid sequence
        variant_type: :obj:`str`
            Type of protein variant. Supported values currently include
            ``"substitution"``, ``"deletion"``, ``"insertion"``, and ``"delins"``.
        sources : :obj:`list` of :obj:`cancermuts.datasources.Datasource`
            sources the variant was derived from
        metadata : :obj:`dict`
            Dictionary encoding variant-associated metadata.

        """
        self.sequence = sequence
        self.start = start
        self.end = end
        self.ref = ref
        self.alt = alt
        self.variant_type = variant_type
        if sources is None:
            self.sources = []
        else:
            self.sources = sources
        if metadata is None:
            self.metadata = {}
        else:
            self.metadata = metadata

    def __str__(self):
        if self.variant_type == "substitution":
            return "%s%d%s" % (self.ref, self.start, self.alt)
        elif self.variant_type == "deletion":
            if self.start == self.end:
                return "%s%ddel" % (self.ref, self.start)
            return "%s%d_%ddel" % (self.ref, self.start, self.end)
        elif self.variant_type == "insertion":
            return "%d_%dins%s" % (self.start, self.end, self.alt)
        elif self.variant_type == "delins":
            return "%s%d_%ddelins%s" % (self.ref, self.start, self.end, self.alt)
        else:
            raise TypeError("Unsupported variant type %s" % self.variant_type)

    def __repr__(self):
        return "<%s %s from %s>" % (self.__class__.__name__, str(self),
                                    ", ".join([s.name for s in self.sources]))         
    def __eq__(self, other):
        if not isinstance(other, ProteinVariant):
            return False
        return (self.sequence == other.sequence and
                self.start == other.start and
                self.end == other.end and
                self.ref == other.ref and
                self.alt == other.alt and
                self.variant_type == other.variant_type)
