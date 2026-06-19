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
import re
import pandas as pd
from Bio.SeqUtils import seq3

class Sequence(object):
    """The most fundamental class of Cancermuts, this class starts from a
    protein sequence definition. It acts as a collection of protein variant objects
    and properties.

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
        Source from where the protein sequence is downloaded
    gene_id : :obj:`str`, optional
        gene name to which the sequence to be downloaded belongs
    sequence : :obj:`str`
        protein sequence as obtained from the data source
    sequence_numbering : :obj:`list of int`
        list of integers, starting from one, each representing a position
    properties : :obj:`dict`
        Dictionary mapping property names to one SequenceProperty object per
        property type. Each SequenceProperty stores its own entries internally
        by residue position
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
        self.variants = VariantRegister()
        self.source = source
        self.sequence = sequence
        self.sequence_numbering = list(range(1, len(self.sequence) + 1))
        self.properties = {}

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
        return iter(zip(self.sequence_numbering, self.sequence))

    def __repr__(self):
        return f"<Sequence gene_id={self.gene_id}, uniprot_ac={self.uniprot_ac}, isoform={self.isoform}, is_canonical={self.is_canonical}, source={self.source.name}, {len(self.sequence)} positions>"

    def add_property(self, property_name, positions, sources=None, **metadata):
        """
        Add an entry to a sequence property.
        A Sequence stores one SequenceProperty object per property name. This method
        creates the property object if needed, then stores the provided entry inside
        that object at the requested residue positions.

        """
        if property_name not in sequence_properties_classes:
            raise ValueError(f"Unknown property_name: {property_name}")

        if positions is None:
            positions = []
        if not isinstance(positions, list):
            raise TypeError("positions must be a list of 1-based residue numbers")
        if len(set(positions)) != len(positions):
            raise ValueError(f"Repeated positions provided for property {property_name}: {positions}")
        if not set(positions).issubset(set(self.sequence_numbering)):
            raise ValueError(f"Property {property_name} is outside sequence bounds for sequence of length {len(self.sequence)}")
        if property_name not in self.properties:
            self.properties[property_name] = sequence_properties_classes[property_name]()

        self.properties[property_name].add_entries(positions=positions, sources=sources, **metadata)

    def _variant_wt_check(self, var):
        if var.start < 1 or var.end > len(self.sequence):
            raise ValueError(f"Variant {var.hgvs} is outside sequence bounds: "
                             f"{var.start}-{var.end} for sequence of length {len(self.sequence)}")
        if self.sequence[var.start - 1:var.end] != var.ref:
            raise ValueError(f"WT mismatch for {var.hgvs}: expected {var.ref} at "
                             f"{var.start}-{var.end}, got {self.sequence[var.start - 1:var.end]}")
        if var.variant_type == "insertion":
            right_pos = var.start + 1
            if right_pos > len(self.sequence):
                raise ValueError(f"Variant {var.hgvs} is outside sequence bounds")
            if self.sequence[right_pos - 1] != var.alt[-1]:
                raise ValueError(f"WT mismatch for {var.hgvs}: expected right flank {var.alt[-1]} at "
                                 f"position {right_pos}, got {self.sequence[right_pos - 1]}")

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
        self._variant_wt_check(var)
        if var.hgvs not in self.variants:
            self.variants.add_variant(var)
            self.log.info("Adding variant %s to sequence %s" % (str(var), self.__repr__()))
        else:
            self.log.info("Variant %s already in sequence %s; will just add sources and metadata" % (str(var), self.__repr__()))
            existing = self.variants[var.hgvs]
            existing.sources.extend(var.sources)

            for k in var.metadata:
                if k in existing.metadata:
                    existing.metadata[k].extend(var.metadata[k])
                    self.log.debug("    metadata %s was extended" % k)
                else:
                    existing.metadata[k] = var.metadata[k]
                    self.log.debug("    metadata %s was added anew" % k)

    def variants_at_position(self, position, variant_types=None):

        if position not in self.sequence_numbering:
            raise ValueError(f"Position {position} is outside sequence bounds: 1-{len(self.sequence)}")
        if variant_types is None:
            variant_types = ProteinVariant.supported_variant_types
        elif isinstance(variant_types, str):
            variant_types = {variant_types}
        else:
            variant_types = set(variant_types)

        invalid_variant_types = variant_types - ProteinVariant.supported_variant_types
        if invalid_variant_types:
            raise ValueError(f"Invalid variant type(s): {invalid_variant_types}")
        matching_variants = []
        for variant in self.variants:
            if variant.variant_type not in variant_types:
                continue
            if variant.variant_type == "insertion":
                if position in {variant.start, variant.start + 1}:
                    matching_variants.append(variant)
            elif variant.start <= position <= variant.end:
                matching_variants.append(variant)
        return matching_variants

    def properties_at_position(self, position, property_names=None):
        """
        Return property entries annotated at a residue position.
        If property_names is provided, only entries for that property type are
        returned. If property_names is None, entries for all property types present
        at the position are returned.
        """
        if position not in self.sequence_numbering:
            raise ValueError(f"Position {position} is outside sequence bounds: 1-{len(self.sequence)}")

        if property_names is None:
            property_names = list(self.properties.keys())
        elif isinstance(property_names, str):
            property_names = [property_names]
        else:
            property_names = list(property_names)
        unknown_properties = set(property_names) - set(sequence_properties_classes)
        if unknown_properties:
                raise ValueError(f"Unknown property name(s): {unknown_properties}")

        properties = {}
        for property_name in property_names:
            if property_name not in self.properties:
                continue
            entries = self.properties[property_name].get_entries_at(position)
            if entries:
                properties[property_name] = entries
        return properties


class ProteinVariant(object):
    """This class describes in-frame protein variants on a reference sequence.

        Attributes
        ----------
        start : :obj:`int`
            1-based start coordinate on the protein sequence
        end : :obj:`int`
            1-based end coordinate on the protein sequence
        ref : :obj:`str`
            reference amino-acid sequence
        alt : :obj:`str`
            altered amino-acid sequence
        sources : :obj:`list` of :obj:`cancermuts.datasources.Datasource`
            sources the variant was derived from
        metadata : :obj:`dict`
            Dictionary encoding variant-associated metadata.
    """
    supported_variant_types = {"missense", "deletion", "insertion", "delins"}
    @logger_init
    def __init__(self, start, end, ref, alt, sources=None, metadata=None):
        """Constructor for the ProteinVariant class.
        
        Parameters
        ----------
        start : :obj:`int`
            1-based start coordinate on the protein sequence
        end : :obj:`int`
            1-based end coordinate on the protein sequence
        ref : :obj:`str`
            reference amino-acid sequence
        alt : :obj:`str`
            altered amino-acid sequence
        sources : :obj:`list` of :obj:`cancermuts.datasources.Datasource`
            sources the variant was derived from
        metadata : :obj:`dict`
            Dictionary encoding variant-associated metadata.

        """
        self.start = start
        self.end = end
        self.ref = ref
        self.alt = alt
        if sources is None:
            self.sources = []
        else:
            self.sources = sources
        if metadata is None:
            self.metadata = {}
        else:
            self.metadata = metadata

        self.variant_type

    @property
    def variant_type(self):
        if (self.start == self.end and len(self.ref) == len(self.alt) == 1
        and self.ref != self.alt):
            return "missense"
        elif self.ref != "" and self.alt == "":
            return "deletion"
        elif (self.start == self.end
        and len(self.alt) > 2
        and self.ref == self.alt[0]):
            return "insertion"
        elif self.ref != "" and self.alt != "" and self.ref != self.alt:
            return "delins"
        else:
            raise ValueError(f"Cannot determine variant type from start={self.start}, "
                             f"end={self.end}, ref='{self.ref}', alt='{self.alt}'")

    @property
    def hgvs(self):
        if self.variant_type == "missense":
            return "p.%s%d%s" % (seq3(self.ref), self.start, seq3(self.alt))
        elif self.variant_type == "deletion":
            if self.start == self.end:
                return "p.%s%ddel" % (seq3(self.ref), self.start)
            return "p.%s%d_%s%ddel" % (seq3(self.ref[0]), self.start, seq3(self.ref[-1]), self.end)
        elif self.variant_type == "insertion":
            return "p.%s%d_%s%dins%s" % (seq3(self.ref), self.start, seq3(self.alt[-1]), self.start + 1, seq3(self.alt[1:-1]))
        elif self.variant_type == "delins":
            if self.start == self.end:
                return "p.%s%ddelins%s" % (seq3(self.ref), self.start, seq3(self.alt))
            return "p.%s%d_%s%ddelins%s" % (seq3(self.ref[0]), self.start, seq3(self.ref[-1]), self.end, seq3(self.alt))

    def __str__(self):
        return self.hgvs

    def __repr__(self):
        return "<%s %s from %s>" % (self.__class__.__name__, str(self),
                                    ", ".join([s.name for s in self.sources]))         
    def __eq__(self, other):
        if not isinstance(other, ProteinVariant):
            return False
        return (self.start == other.start and
                self.end == other.end and
                self.ref == other.ref and
                self.alt == other.alt)

class VariantRegister:
    """
    Thin abstraction layer over a dictionary storing ProteinVariant objects
    with their HGVS string representation as key.
    """
    _aa1 = r"[ACDEFGHIKLMNPQRSTVWY]"
    _aa3 = r"(?:Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Ter)"

    _mavisp_formats = [rf"{_aa1}\d+{_aa1}",
                       rf"{_aa1}\d+del",
                       rf"{_aa1}\d+_{_aa1}\d+del",
                       rf"{_aa1}\d+_{_aa1}\d+ins{_aa1}+",
                       rf"{_aa1}\d+delins{_aa1}+",
                       rf"{_aa1}\d+_{_aa1}\d+delins{_aa1}+"]

    _hgvs_formats = [rf"p\.{_aa3}\d+{_aa3}",
                     rf"p\.{_aa3}\d+del",
                     rf"p\.{_aa3}\d+_{_aa3}\d+del",
                     rf"p\.{_aa3}\d+_{_aa3}\d+ins(?:{_aa3})+",
                     rf"p\.{_aa3}\d+delins(?:{_aa3})+",
                     rf"p\.{_aa3}\d+_{_aa3}\d+delins(?:{_aa3})+"]

    _mavisp_format = re.compile(rf"({'|'.join(_mavisp_formats)})")
    _hgvs_format = re.compile(rf"({'|'.join(_hgvs_formats)})")

    def __init__(self):
        self._register = {}

    def add_variant(self, variant):
        """
        Add a ProteinVariant to the register.

        Parameters
        ----------
        variant : :obj:`cancermuts.core.ProteinVariant`
            Variant to be added to the register.
        """
        if not isinstance(variant, ProteinVariant):
            raise TypeError("VariantRegister only accepts ProteinVariant objects, "
                           f"got {type(variant)}")
        hgvs = variant.hgvs
        if hgvs in self._register:
            raise ValueError(f"Variant {hgvs} already present in register")
        self._register[hgvs] = variant

    def _is_valid_mavisp_string(self, item):
        return self._mavisp_format.fullmatch(item) is not None

    def _is_valid_hgvs_string(self, item):
        return self._hgvs_format.fullmatch(item) is not None
    
    def _mavisp2hgvs(self, item):
        item = re.sub(r"([A-Z])(\d+)", lambda m: f"{seq3(m.group(1))}{m.group(2)}", item)
        item = re.sub(r"(ins|delins)([A-Z]+)$", lambda m: f"{m.group(1)}{seq3(m.group(2))}", item)
        item = re.sub(r"(\d+)([A-Z])$", lambda m: f"{m.group(1)}{seq3(m.group(2))}", item)
        return f"p.{item}"

    def __getitem__(self, item):
        """
        Retrieve a variant from the register by one-letter or HGVS string.

        Parameters
        ----------
        item : :obj:`str`
            Variant string in MAVISp one-letter format (e.g. ``H39F``) or HGVS protein format
            (e.g. ``p.His39Phe``).

        Returns
        -------
        :obj:`cancermuts.core.ProteinVariant`
            The corresponding stored variant.
        """
        if self._is_valid_hgvs_string(item):
            return self._register[item]
        elif self._is_valid_mavisp_string(item):
            return self._register[self._mavisp2hgvs(item)]
        raise KeyError(f"Unsupported variant format: {item}")

    def _sorted_variants(self):
        return sorted(self._register.values(),
        key=lambda var: (var.start, var.end, var.variant_type, var.ref, var.alt))

    def __iter__(self):
        """
        Iterate over ordered stored variants.
        """
        return iter(self._sorted_variants())

    def __contains__(self, item):
        try:
            self[item]
            return True
        except KeyError:
            return False

    def get_variant_table(self, variant_types=None, metadata=[], metadata_formatter=None):
        """
        Return a pandas DataFrame containing the ordered variants stored in the register.

        Parameters
        ----------
        variant_types : :obj:`str` or iterable of :obj:`str`, optional
            Variant type to include in the output. If provided, only variants whose
            ``variant_type`` matches one of the specified values are returned.
            If None, all variants are included.
        metadata : iterable of :obj:`str`, optional
            Metadata fields to include as additional columns.
        metadata_formatter : callable, optional
            Function used to format metadata values before adding them to the table.

        Returns
        -------
        :obj:`pandas.DataFrame`
            DataFrame with one row per variant
        """

        variants = self._sorted_variants()
        if variant_types is not None:
            if isinstance(variant_types, str):
                variant_types = {variant_types}
            else:
                variant_types = set(variant_types)
            invalid_var_types = variant_types - ProteinVariant.supported_variant_types
            if invalid_var_types:
                raise ValueError(f"Invalid variant type(s): {invalid_var_types}")
            variants = [v for v in variants if v.variant_type in variant_types]

        variant_dictionary = {"variant_hgvs": [],
                              "variant_start": [],
                              "variant_end": [],
                              "variant_ref": [],
                              "variant_alt": [],
                              "variant_type": [],
                              "sources": []}

        for md in metadata:
            variant_dictionary[md] = []

        for v in variants:
            variant_dictionary["variant_hgvs"].append(v.hgvs)
            variant_dictionary["variant_start"].append(v.start)
            variant_dictionary["variant_end"].append(v.end)
            variant_dictionary["variant_ref"].append(v.ref)
            variant_dictionary["variant_alt"].append(v.alt)
            variant_dictionary["variant_type"].append(v.variant_type)

            source_names = sorted({source.name for source in v.sources})
            variant_dictionary["sources"].append(",".join(source_names) if source_names else None)
            for md in metadata:
                if md in v.metadata:
                    values = v.metadata[md]
                else:
                    values = None
                if metadata_formatter is not None:
                    values = metadata_formatter(values, md)
                variant_dictionary[md].append(values)

        return pd.DataFrame(variant_dictionary)
