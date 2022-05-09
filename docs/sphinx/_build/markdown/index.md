<!-- cancermuts documentation master file, created by
sphinx-quickstart on Mon May  2 10:40:47 2022.
You can adapt this file completely to your liking, but it should at least
contain the root `toctree` directive. -->
# API Reference


* [Index](genindex.md)


* [Module Index](py-modindex.md)

## Core cancermuts classes — `cancermuts.core`

core classes that define the basic framework to handle a sequence
and its mutations.


### _class_ cancermuts.core.Mutation(sequence_position, mutated_residue_type, sources=None, metadata=None)
This class describe a missense mutation as amino-acid replacement.


#### sequence_position()
SequencePosition object to which this mutation belongs, corresponding
to the residue that is mutated


* **Type**

    `cancermuts.core.SequencePosition`



#### sources()
source the mutation was derived from


* **Type**

    `list` of `cancermuts.datasources.Datasource`



#### mutated_residue_type()
single-letter code for the mutated residue type for this position


* **Type**

    `str`



#### mutations()
list of mutations for this sequence position (if any)


* **Type**

    `list` of cancermuts.core.Mutation objects



#### properties()
dictionary encoding position properties. This dictionary needs to have
`str` as key and cancermuts.properties.PositionProperty as value.


* **Type**

    `dict`



#### \__init__(sequence_position, mutated_residue_type, sources=None, metadata=None)
Constructor for the Mutation class.


* **Parameters**

    
    * **sequence_position** (`cancermuts.core.SequencePosition`) – SequencePosition object to which this mutation belongs, corresponding
    to the residue that is mutated
    single-letter code wild-type residue for this position


    * **mutated_residue_type** (`str`) – single-letter code for the mutated residue type for this position


    * **sources** (`list` of `cancermuts.datasources.Datasource`) – list of one or more sources the mutation was derived from


    * **metadata** (`dict` or `None`) – dictionary encoding position properties. This dictionary needs to have
    `str` as key and cancermuts.properties.PositionProperty as value.



### _class_ cancermuts.core.Sequence(gene_id, sequence, source, aliases=None)
The most fundamental class of Cancermuts, this class starts from a
protein sequence definition. It acts as a collection of ordered
SequencePosition objects which depend on the specific
sequence itself.


* **Parameters**

    
    * **gene_id** (`str`) – ID of the gene (gene name) to which the sequence belongs


    * **sequence** (str) – Protein sequence corresponding to the gene and isoform of interest


    * **source** (`cancermuts.datasources.DataSource`) – source from where the Sequence has been downloaded from


    * **aliases** (`dict` of `str`, both for key and value, or None, optional) – This argument assigns the aliases dictionary for the Sequence class.
    if None, an empty dictionary is created.



#### source()
Source from where the protein sequence is downloaded from


* **Type**

    `cancermuts.datasources.Datasource`



#### gene_id()
gene name to which the sequence to be downloaded belongs to


* **Type**

    `int`, optional



#### sequence()
protein sequence as obtained from the data source


* **Type**

    `str`



#### sequence_numbering()
list of integers, starting from one, each representing a position


* **Type**

    `list of int`



#### properties()
dictionary including the downloaded protein-associated properties.
These can span one or more residues.


* **Type**

    `dict`



#### aliases()
This dictionary maps different “aliases” (i.e. identifiers) for the
protein of interest. The key should represent the type of identifier
it’s being stored in the value and can be any string. For instance,
one could have “cosmic” : “AMBRA1” to mean that the protein AMBRA1
has identifier AMBRA1 in cosmic


* **Type**

    `dict`



#### add_property(prop)
Adds sequence property to sequence object. If a property of the same
category is already present, the property will be added to the same
category; otherwise the category will be created anew


* **Parameters**

    **prop** (`cancermuts.properties.SequenceProperty`) – new property to be added



#### index2seq(idx)
Returns the sequence numbering corresponding to the idx-th
residue, starting from 0. usually this corresponds to idx+1.


* **Parameters**

    **idx** (`int`) – 0-index position in the protein sequence



* **Returns**

    sequence number corresponding to the idx-th position



* **Return type**

    i:obj:int



#### seq2index(seqn)
Returns the 0-indexed position number in the sequence corresponding
to sequence number seqn. It is usually equal to seqn-1..


* **Parameters**

    **seqn** (`int`) – residue position number in the protein sequence



* **Returns**

    sequential number corresponding to the position in the sequence
    (starting from 0)



* **Return type**

    `int`



### _class_ cancermuts.core.SequencePosition(wt_residue_type, sequence_position, mutations=None, properties=None)
This class describe a certain sequence position in a protein.
SequencePositions usually belong to a Sequence object.
It is possible to annotate a Sequence position with either a Mutation
or a PositionProperty.


#### source()

* **Type**

    `cancermuts.datasources.Datasource`



#### wt_residue_type()
single-letter code wild-type residue for this position


* **Type**

    `str`



#### sequence_position()
number corresponding to the sequence position for this position


* **Type**

    `int`



#### mutations()
list of mutations for this sequence position (if any)


* **Type**

    `list` of cancermuts.core.Mutation objects



#### properties()
dictionary encoding position properties. This dictionary needs to have
`str` as key and cancermuts.properties.PositionProperty as value.


* **Type**

    `dict`



#### \__init__(wt_residue_type, sequence_position, mutations=None, properties=None)
Constructor for the SequencePosition class.


* **Parameters**

    
    * **wt_residue_type** (`str`) – single-letter code wild-type residue for this position


    * **sequence_position** (`int`) – number corresponding to the sequence position for this position


    * **mutations** (`list` of cancermuts.core.Mutation objects or None) – list of mutations for this sequence position to be added


    * **properties** (`dict` or `None`) – dictionary encoding position properties. This dictionary needs to have
    `str` as key and cancermuts.properties.PositionProperty as value.



#### add_mutation(mut)
Adds mutation to a `SequencePosition` object. If the mutation (in
terms of amino-acid substitution) is already present, source and
metadata will be added to the already-present mutation. Otherwise
the mutation is added anew.


* **Parameters**

    **mut** (`cancermuts.core.Mutation`) – new Mutation object to be added to the SequencePosition



#### add_property(prop)
Adds position property to a `SequencePosition` object. If a
property of the same category is already present, it will be overwritten.
Otherwise, the property is just added to the object.


* **Parameters**

    **prop** (`cancermuts.properties.PositionProperty`) – new Mutation object to be added to the SequencePosition


## datasources classes — `cancermuts.datasources`

Classes to interrogate data sources and annotate various data


### _class_ cancermuts.datasources.COSMIC(database_files=None, database_encoding=None, cancer_type=None)
Bases: `cancermuts.datasources.DynamicSource`, `object`


### _class_ cancermuts.datasources.DynamicSource(\*args, \*\*kwargs)
Bases: `cancermuts.datasources.Source`, `object`

Base class for implementing dynamic data sources. Dynamic data sources
are remote data sources that can be queried through the internet.


### _class_ cancermuts.datasources.ELMDatabase()
Bases: `cancermuts.datasources.DynamicSource`, `object`


### _class_ cancermuts.datasources.ELMPredictions()
Bases: `cancermuts.datasources.DynamicSource`, `object`


### _class_ cancermuts.datasources.ManualAnnotation(datafile, \*\*parsing_options)
Bases: `cancermuts.datasources.StaticSource`


### _class_ cancermuts.datasources.MobiDB()
Bases: `cancermuts.datasources.DynamicSource`


### _class_ cancermuts.datasources.MyVariant()
Bases: `cancermuts.datasources.DynamicSource`, `object`


### _class_ cancermuts.datasources.PhosphoSite(database_files=None)
Bases: `cancermuts.datasources.DynamicSource`, `object`


### _class_ cancermuts.datasources.Source(name, version, description)
Bases: `object`

Base class for implementing data sources.


#### name()
name of the data source


* **Type**

    `str`



#### version()
version of the data source


* **Type**

    `str`



#### description()
short description of the data source


* **Type**

    `str`



#### \__init__(name, version, description)
Class constructor


* **Parameters**

    
    * **name** (`str`) – name of the data source


    * **version** (`str`) – version of the data source


    * **description** (`str`) – short description of the data source



### _class_ cancermuts.datasources.StaticSource(\*args, \*\*kwargs)
Bases: `cancermuts.datasources.Source`, `object`

Base class for implementing static data sources. Static data sources
are local to the system in use and need to be provided manually. They
don’t change without manual intervention


### _class_ cancermuts.datasources.UniProt(\*args, \*\*kwargs)
Bases: `cancermuts.datasources.DynamicSource`, `object`

Class for the UniProt data source. It is used to download the protein
sequence of the main UniProt isoform for a certain gene and build the
Sequence object, which is the main entry point for annotations in cancermuts


### _class_ cancermuts.datasources.cBioPortal(cancer_studies=None, max_connections=10000, max_threads=40)
Bases: `cancermuts.datasources.DynamicSource`, `object`


### _class_ cancermuts.datasources.gnomAD(version='2.1')
Bases: `cancermuts.datasources.DynamicSource`, `object`

## metadata classes — `cancermuts.metadata`

Classes to handle metadata

## datasources classes — `cancermuts.table`

Classes to save/read the table and generate a graphical representation

## metadata classes — `cancermuts.log`

Logging facilities for cancermuts
