# Tutorial

This tutorial will show an example showcasing the usage of Cancermuts on a
well-known autophagy marker protein, MAP1LC3B (or more simply, LC3B).

We will download associated cancer mutations from both cBioPortal and COSMIC
and annotate proteins and mutations with all the data sources available in
Cancermuts. 

## Preliminary operations

This tutorial expects that the Cancermuts Python package has been installed and
is available to be imported. If you have followed the recommended installation
procedure and installed Cancermuts in a Python environment you need to activate
it first, as detailed in point 3 of the installation.

You also need to have downloaded the COSMIC mutation export file, as detailed in 
the installation section.

Finally, you can decide to perform this tutorial either interactively (i.e. 
running line by line in a Python command line) or as a script (i.e. saving
the lines you need in a script and running it). The two options are equvalent
from the point of view of the result. If you're going for the first
option, we recommended that you install the IPython command line interface
which is more user friendly than the regular Python interpreter:

```
$ pip install ipython
```

and then just run `ipython` to enter the interpreter:

```(cancermuts_env) teo@kb-bioinfo01:cancermuts$ ipython
Python 3.8.0 (default, Dec  9 2021, 17:53:27)
Type 'copyright', 'credits' or 'license' for more information
IPython 8.3.0 -- An enhanced Interactive Python. Type '?' for help.

In [1]:
```

lastly, you can check if your Cancermuts installation was successful by
trying and import the cancermuts package:

```
In [1]: import cancermuts

In [2]:
```

you shouldn't receive any error. If you receive an error instead we recommend
checking that your installation was performed correctly and/or open an Issue
on our GitHub repository to receive assistance.

## Introduction

Cancermuts works by interrogating different *data sources* to gather the
requested information. These are all encoded in the `datasources` module
which includes one class for every data source available in cancermuts.
The steps for obtaining a cancermuts annotation are roughly as follows:

1. Create a Sequence object by interrogating an appropriate data source.
The Sequence object is the main container of information in cancermuts
and links together all the available annotations

2. Download mutations and add them to the sequence object. This is done by
appropriately using other data sources

3. Download annotations for either mutations (these are called metadata) or
position or sequence properties (these are called properties). This is once
again done by using appropriate data sources for the task.

4. Create a pandas data frame from the Sequence object. This will contain all
the gathered information.

5. Create a visualization from the data frame

## Tutorial script

### Canonical isoform tutorial

The steps performed in the following tutorial are also written in a pre-built
Python script available in the `docs` folder, called `tutorial.py`. In order
to run it, you should have installed the Cancermuts package and activated the
virtual environment in which it is installed (see Installation). This tutorial
script follows the tutorial steps and has been written using a single cancer
type as reference (colorectal cancer).

{% hint style='tip' %}
Depending on the location of the files required by Cancermuts on your system,
you will need to edit the hardcoded paths in the script to match what is
available to you
{% endint %}

Once that is done, from the `cancermuts/docs` directory, you can just run:

```
$ python tutorial.py
```

The result should be a `metatable.csv` output file together with a `my_table.pdf`
figure.

We also provide a similar example which considers all available cancer types,
in the `tutorial_pancancer.py` file. It's run in the same way as with the
other tutorial script:

```
$ python tutorial_pancancer.py
```

Similarly as before, the result should be a `metatable_pancancer.csv` output file
together with a `my_table_pancancer.pdf` figure.

### Isoform tutorial

Cancermuts allows requesting alternative isoforms, if available in UniProt.
To request a specific isoform, provide a isoform ID, in the form of a UniProt isoform
identifier, to the isoform argument, as in the following example. If this is not done,
the canonical UniProt isoform will be used. Notice that, currently, only some
data sources support alternative isoform; those that do not support them will
raise exceptions if a Sequence object containing an alternative isoform is provided to them.
In this case we provide an example which use a non-canonical isoform as input in the 
`tutorial_isoforms.py` file. It's run in the same way as with the other 
tutorial scripts:

```
$ python tutorial_isoforms.py
```
It loads a specific AMBRA1 isoform (here we use AMBRA1 because LC3B has not-main isoforms) and 
ends after successfully downloading and displaying the isoform sequence.

## Tutorial steps

### The Sequence object

Cancermuts works 
The first operation to be performed when starting to work with Cancermuts is creating
a Sequence object. This object represents the protein sequence we want to consider for
the analysis and links together all the information we are going to collect. In order to
create it we use a data source class for UniProt from which we can downlaod protein sequences.
This is only source for sequences available at the moment, more can be added in the future.

To download the canonical UniProt sequence, you will need the corresponding HGCN symbol only.
The UniProt object automatically identifies the most likely UniProt ID and AC corresponding to the
gene name. Alternatively, it is possible to provide the UniProt ID manually.

This is done as follows:


```py
# import the UniProt data source class
>>> from cancermuts.datasources import UniProt

# create the corresponding uniprot object
>>> up = UniProt()

# get the sequence for the protein
>>> seq = up.get_sequence('MAP1LC3B')

# alternatively, we can specifically ask for a Uniprot ID and/or AC
>>> seq = up.get_sequence('MAP1LC3B', upid='MLP3B_HUMAN', upac='Q9GZQ8')

# this prints the downloaded protein sequence
>>> print(seq.sequence)
MPSEKTFKQRRTFEQRVEDVRLIREQHPTKIPVIIERYKGEKQLPVLDKTKFLVPDHVNMSELIKIIRRRLQLNANQAFFLLVNGHSMVSVSTPISEVYESEKDEDGFLYMVYASQETFGMKLSV

# the seq.positions attribute is an ordered list of the protein positions:
>>> seq.positions[0:3]

```

```py

# OPTIONAL: Load a specific UniProt isoform instead of canonical. In this example
# we use the AMBRA1 protein:
>>> isoform_id = "Q9C0C7-3"
>>> up = UniProt()
>>> seq = up.get_sequence("AMBRA1", isoform=isoform_id)
>>> print(seq.isoform)
Q9C0C7-3
>>> print(seq.is_canonical)
False
```

### Collecting cancer mutations

Now that we have the Sequence ready, we can start adding annotations to it.
We will first download cancer mutations from cBioPortal and COSMIC.

We first import the required datasource classes for them:

```py
# import data sources classes
>>> from cancermuts.datasources import cBioPortal, COSMIC
```

#### cBioPortal

For cBioPortal, we first create the respective source object:

```py
>>> cb = cBioPortal(cancer_studies=['coadread_dfci_2016', 
	                            'coadread_genentech',
	                            'coadread_tcga_pan_can_atlas_2018'])
```

For this tutorial, we specify a list of cancer studies Cancermuts should 
use, which correspond to a few colorectal cancer studies. This is especially
useful if we are interested in e.g. a certain cancer type or some other
specific studies. Studies can be referred to using their cBioPortal 
study identifier.

{% hint style='info' %}
Finding the correct set of identifier for the studies of interest can
be tricky, as it isn't obvious from the cBioPortal website which study
IDs are connected to specific studies. When creating a cBioPortal object
without the `cancer_studies` argument, two attributes of the object can
be helpful to this regard. `cb.cancer_types` is a list of cancer types
as classified by cBioPortal and `cb.cancer_studies` is the full list of
studies available on cBioPortal. These list are regular pandas dataframes
and can be manipulated as such, including saved to csv files or similar.

Alternatively, one can visit the [cBioPortal datasets page](https://www.cbioportal.org/datasets)
which contains a list of all the available studies. While this list doesn't
contain a column with study IDs, the web links that link the studies in the
"Name" column do contain the IDs. For instance, the study named 
"Breast Invasive Carcinoma (Broad, Nature 2012)" links the page
`https://www.cbioportal.org/study?id=brca_broad` - meaning the corresponding
study ID is `brca_broad`.
{% endhint %}

{% hint style='info' %}
It is also possible to initialize the cBioPortal data source object without
specifying cancer studies:

```py
>>> cb = cBioPortal()
```

In this case, Cancermuts will gather mutations from all cancer studies available
in cBioPortal.
{% endhint %}

Finally, we use this object to gather the mutations from the selected cancer
studies and add them to our Sequence object we created earlier:

```py
>>> cb.add_mutations(seq)
```

It is possible to downloaded metadata about the downloaded aminoacid mutations
as well by using the `metadata` argument of `add_mutations`, which supports
a list of strings, one for each metadata type. This is usually recommended.
By default no metadata are added. Supported metadata for cBioPortal are:

* `cancer_type`: type of cancer the mutation was found in, depending on the study
* `cancer_study`: cBioPortal study the mutation was found in
* `genomic_coordinates`: Genomic coordinates of the corresponding genomic mutation
* `genomic_mutations`: Genomic coordinates and base pair substitution
    of the corresponding genomic mutation

So for instance:

```py
>>> cb.add_mutations(seq, metadata=['cancer_type', 'cancer_study', 'genomic_mutations'])
```

{% hint style='info' %}
the gene or protein name is not among the arguments - this is
because it is inferred from the `seq` object. Only mutations
of the corresponding gene will be annotated.
{% endhint %}

This step might take a few minutes as Cancermuts interrogates the cBioPortal
online database. The final result of this step will have modified the `seq`
Sequence object we have created by adding mutations to the positions of the
respective residues. In this case, Cancermuts identified only three
mutations for cBioPortal:

```py
>>> print(seq.positions[38].mutations)
[<Mutation K39R from cBioPortal>]

>>> print(seq.positions[64].mutations)
[<Mutation K65E from cBioPortal>]

>>> print(seq.positions[122].mutations)
[<Mutation L123S from cBioPortal>]
```

Each mutation is recorded in a Mutation object that can be further explored:

```py
>>> seq.positions[64].mutations
[<Mutation K65E from cBioPortal>]

>>> seq.positions[64].mutations[0].sources
[<cancermuts.datasources.cBioPortal at 0x1510a0455c50>]

>>> seq.positions[64].mutations[0].mutated_residue_type
'E'

>>> seq.positions[38].mutations[0].metadata
{'cancer_type': [<CancerType Colorectal Adenocarcinoma from cBioPortal>],
 'cancer_study': [<CancerStudy coadread_genentech from cBioPortal>],
 'genomic_mutations': [<GenomicMutation hg19,16:g.87435877A>G from cBioPortal>],
 'genomic_coordinates': [<GenomicCoordinates 16:87435877-87435877 in hg19 from cBioPortal>]}
 ```

{% hint style='info' %}
The cancer type is not always present in the data downloaded from cBioPortal.
In this cases, Cancermuts tries its best to infer it from the study name.
{% endhint %}

#### COSMIC

As before, we first create a COSMIC data source object:

```py
cosmic = COSMIC(database_files='/data/databases/cosmic-v95/CosmicMutantExport.tsv',
	            database_encoding=['latin1'])
```

here the `database_files` argument is a string. If the user wants to use
more than one database file, these should be provided as a list of strings.
Usually, the argument for this file would be the COSMIC database file  that 
was downloaded as detailed in the Install section.

Similarly, `database_encoding` defines the
text file encoding for every file (it is `latin1` for COSMIC version 95).

{% hint style='danger' %}
As the default database file is rather large, we recommend running this step on
a computer with at least 32 GB of free memory. Otherwise, it is possible to
filter the database file first, for instnace keeping only the rows that
contain the gene name of interest. In bash this can be done by running:

```
$ head -n 1 CosmicMutantExport.tsv > header.txt
$ grep MAP1LC3B CosmicMutantExport.tsv > content.txt
$ cat header.txt content.txt > COSMIC_map1lc3b.csv
$ rm header.txt content.txt
```

and then using the resulting file as the database file
{% endhint %}

We can then add mutations from the COSMIC data source:

```py
cosmic.add_mutations(seq, 
					 cancer_sites=['large_intestine'],
					 cancer_site_subtype_1=['colon'],
					 cancer_types=['carcinoma'],
					 cancer_histology_subtype_1=['adenocarcinoma'], 
					 metadata=['genomic_coordinates', 'genomic_mutations',
					 		   'cancer_site', 'cancer_histology'])
```

Here we restrict the search to those mutations that are involved in colorectal
adenocarcinoma. 

{% hint style='info' %}
It is also possible to search in any available cancer type or site as well (i.e. pancancer),
by not specifying cancer types or sites in the `add_mutations` call, for instance:

```py
cosmic.add_mutations(seq, metadata=['genomic_coordinates', 'genomic_mutations', 
                                       'cancer_site', 'cancer_histology'])
```
{% endhint %}

It should be noted that COSMIC supports up to four cancer histology types
and four cancer histology subtypes and Cancermuts allows to filter by any of this.
In particular, only the mutations that  correspond to *all* the selected criteria
are retained. Please see the [COSMIC phenotype classification](https://cancer.sanger.ac.uk/cosmic/classification)
for the classification that COSMIC uses.

Furthermore, similarly to cBioPortal, we retain some metadata:

* `cancer_histology`: histology information on the cancer(s) in which this mutation was found
* `cancer_site`: site information on the cancer(s) in which this mutation was found
* `genomic_coordinates`: Genomic coordinates of the corresponding genomic mutation
* `genomic_mutations`: Genomic coordinates and base pair substitution
    of the corresponding genomic mutation

In this case, filtering the database allows to add only a single mutation,
which is K65E. It should be noted that the same amino acid substitution 
was already identified by cBioPortal. This means that the mutation object
corresponding to this aminoacid substitution is annotated with metadata for
this newfound mutation. We see now that the details obtained by both cBioPortal
and COSMIC are present:

```py
>>> seq.positions[64].mutations[0]
<Mutation K65E from cBioPortal,COSMIC>

>>> seq.positions[64].mutations[0].sources
[<cancermuts.datasources.cBioPortal at 0x1510a0455c50>,
 <cancermuts.datasources.COSMIC at 0x15109ff10050>]

>>> seq.positions[64].mutations[0].metadata
{'cancer_type': [<CancerType Colorectal Adenocarcinoma from cBioPortal>],
 'cancer_study': [<CancerStudy coadread_genentech from cBioPortal>],
 'genomic_mutations': [<GenomicMutation hg19,16:g.87435877A>G from cBioPortal>,
  <GenomicMutation hg38,16:g.87402271A>G from COSMIC>],
 'genomic_coordinates': [<GenomicCoordinates 16:87435877-87435877 in hg19 from cBioPortal>,
  <GenomicCoordinates 16:87402271-87402271 in hg38 from COSMIC>],
 'cancer_site': [<CancerSite, large_intestine, colon>],
 'cancer_histology': [<TumorHistology, carcinoma, adenocarcinoma>]}
```

### Additional mutation metadata

Cancermuts allows to enrich the downloaded mutations with further metadata. We
will see now how to download the REVEL pathogenicy score and the gnomAD allele
frequency for each specific variant.

#### REVEL score from MyVariant

We will download scores for the REVEL predictor of pathogenicity from the
[MyVariant.info](https://myvariant.info) database. REVEL is associated to
genomic mutations and not to amino-acid mutations, therefore our amino acid
mutations need to be annotated with the genomic mutations metadata in order
for this to be possible. This is performed as detailed above. In order to
download the associated REVEL score, we can use the appropriate data source
class:

```py
>>> from cancermuts.datasources import MyVariant

>>> mv = MyVariant()
>>> mv.add_metadata(seq)
```

we can check that the mutations have been annotated with the REVEL score:

```py
>>> seq.positions[64].mutations[0].metadata['revel_score']
[<DbnsfpRevel, 0.314>, <DbnsfpRevel, 0.314>]
```

In this case we see the same score twice, as the mutation was previously
annotated with two genomic mutations. The two genomic mutations corresponded
to the same mutations annotated in two different genome assemblies, therefore
the two scores we are able to gather have the same value.

#### gnomAD allele frequencies

Similarly, we annotate mutations with their exome or genome allele frequencies
as found in the [gnomAD database](https://www.gnomad.org). It is also possible to
annotate with the highest allele frequency found in non-bottlenecked 
populations (referred to as popmax). This works as you would expect by now:

```py
>>> from cancermuts.datasources import gnomAD

>>> gnomad = gnomAD(version='2.1')
>>> gnomad.add_metadata(seq, md_type=['gnomad_exome_allele_frequency',
	                              'gnomad_genome_allele_frequency',
                                  'gnomad_popmax_exome_allele_frequency',
	                              'gnomad_popmax_genome_allele_frequency'])
```

here, we specify the `version` argument to specify the version of gnomAD to be 
considered. Please refer to the API documentation for all the available versions.

The `md_type` keyword allows to select which metadata type(s) to annotate, i.e.
choose between exome or genome allele frequency (or both as in the example) and whether 
to annotate the popmax as well. Note, that they can all be chosen independently 
of each other.

We calculate the allele frequency as the ratio between the total allele count over
the allele number as found in gnomAD, if the entry for the corresponding variant 
is available.

The popmax allele frequency is found by calculating the exome and/or genome allele frequency from 
all supported non-bottlenecked populations and then finding the maximum frequency of those.

Finally we can check the downloaded metadata:

```py
>>> seq.positions[64].mutations[0].metadata['gnomad_exome_allele_frequency']
[<gnomADExomeAlleleFrequency, 0.000028>,
 <gnomADExomeAlleleFrequency, 0.000028>]

>>> seq.positions[64].mutations[0].metadata['gnomad_genome_allele_frequency']
[<gnomADGenomeAlleleFrequency, nan>,
 <gnomADGenomeAlleleFrequency, nan>]
```

here Cancermuts could assign an exome allele frequency for the associated genomic
mutations, however it wasn't able to download or calculate the corresponding
genome allele frequency.

### Position properties

We will further annotate our sequence object with properties connected to each single
residue. For the moment, these are post-translation modifications from Phosphosite and
order/disorder propensity from MobiDB.

We import the relative data source class, similarly as what done previously:

```py
from cancermuts.datasources import PhosphoSite, MobiDB
```

#### Post-translational modifications with phosphosite

We will annotate post-translational modifications identified in experiments from
the PhosphoSite Plus database. In order to do so we need a local copy of the
Phosphosite dataset - please see the instructions in the installation guide.

We first create the Phosphosite data source object. We will need to supply the
location of the database files:

```py
>>> ps = PhosphoSite('/data/databases/phosphosite/')
```

by default, Cancermuts expect the file names in the database to be the default
ones. However it's possible to specify the file name for each post-translational
modifications by supplying the `database_files` argument, which should be a
dictionary associating each PTM to a file name:

```py
>>> my_databse_files = {  'acetylation'     : 'my_Acetylation_site_dataset',
                          'methylation'     : 'my_Methylation_site_dataset',
                           ... }
```

the keys of the dictionaries should be the options supported in the `properties`
argument of the `add_position_properties` (see below).

Once the object is created we can add the position properties to our sequence
object. If no `properties` is supplied, all of them will be considered. Otherwise,
the `properties` object should be a list of the keywords corresponding to the PTMs 
that we want to be annotated, as follows:

| Keyword | Description |
| ------- | ----------- |
| `acetylation` | Acetylation | 
| `methylation` | Methylation |
| `O-GalNAc` | O-linked β-N-acetylglucosamine |
| `O-GlcNAc` | O-linked α-N-acetylgalactosamine |
| `phosphorylation` | Phosphorylation |
| `sumoylation` | Sumoylation |
| `ubiquitination` | Ubiquitination |

So, for instance:

```py
>>> ps.add_position_properties(seq, 
	                           properties=['phosphorylation', 'ubiquitination'])
```

However, in this case we will consider all of them:

```py
>>> ps.add_position_properties(seq)
```

Once again, we can check the result:

```py
>>> seq.positions[4].properties
{'ptm_ubiquitination': <PositionProperty Ubiquitination Site from PhosphoSite>}

>>> seq.positions[28].properties
{'ptm_phosphorylation': <PositionProperty Phosphorylation Site from PhosphoSite>}
```

#### structured regions with MobiDB

We use the MobiDB website to annotate disordered regions in our
protein of interest using the consensus of the entry. This works similarly as before:

```py
mdb = MobiDB()
mdb.add_position_properties(seq)
```

We can then check the annotation as done previously:

```py
>>> seq.positions[0].properties['mobidb_disorder_propensity']
<StructuralDisorder, None>

>>> seq.positions[10].properties['mobidb_disorder_propensity']
<StructuralDisorder, Structured, observed>
```

There are four types of annotations, here listed in order of priority.
* Disordered, curated
    * Based on manually curated data
* Disordered, derived
    * Based on primary data, e.g. PDB structures
* Disordered, homology
    * Based on homology inference
* Disordered, predicted
    * Based on predictions

If none of these data are available for the consensus, the residue will not be annotated. 
For more information about the types of evidence, please see the [MobiDB vocabulary](https://mobidb.bio.unipd.it/about/vocabulary)

### Sequence properties

We further annotate properties to the sequence, i.e. properties that cover
multiple residues. We first import the respective classes:

```py
>>> from cancermuts.datasources import ELMPredictions
```

#### Predictions of short linear motifs from ELM

We annotate the sequence using predictions for short linear motifs from the
[Eukaryotic Linear Motif](http://elm.eu.org) database:

```py
>>> elm = ELMPredictions()
>>> elm.add_sequence_properties(seq,
	                            exclude_elm_classes="MOD_.")
```

here `exclude_elm_classes` is a single regular expression that allows to remove
specific ELM classes by specifying a regular expression that matches one or more
[ELM class identifier](http://elm.eu.org/elms). Here we use `"MOD_."` to exclude
post-translational modifications which are already supplied by PhosphoSite.

We can check the linear motifs we have collected:

```py
>>> seq.properties
{'linear_motif': [<SequenceProperty Linear motif from ELM, positions 10,11,12>,
  <SequenceProperty Linear motif from ELM, positions 69,70,71>,
  <SequenceProperty Linear motif from ELM, positions 5,6,7,8,9>,
  ...
```

and their specifics:

```py
>>> seq.properties['linear_motif'][0].type
'NRD cleavage site'

>>> seq.properties['linear_motif'][0].positions
[<SequencePosition, residue R at position 10>,
 <SequencePosition, residue R at position 11>,
 <SequencePosition, residue T at position 12>]
```

Alternatively, and in a very similar fashion, it is possible to use the [gget Python
package](https://github.com/pachterlab/gget) to obtain short linear motifs definitions.
The current implementation only considers "regexp" type of predictions from the
full protein sequence. It should be noted that, unlike when using the ELM webserver
as detailed above to obtain these data, no filtering is applied. For instance:

```py
>>> elm = ggetELMPredictions()
>>> elm.add_sequence_properties(seq,
	                            exclude_elm_classes="MOD_.")
```

This can be useful when e.g. the ELM webserver is not available or for when
many calls in a row are necessary. The ELM webserver requires a minimum 3-minute
interval between queries which has been baked into the current implementation of
the ELM data source.

### Custom annotations

We can further add annotations manually to our dataset. This is for data that
is not available in the databases as of yet, but is useful to have annotated
in the pipeline to have a complete picture. This can be done by using a custom
CSV file as data source. The CSV can contain different types and levels of
information. The file should have the following columns, separated by `;`:

| Column | Description |
| ------ | ----------- |
| `name` | name that is given to this feature. It can be any string. | 
| `site` | position of this feature; see below | 
| `type` | type of this feature; see below |
| `function` | function description of this feature or functional annotation; see below |
| `reference` | reference to the literature for this feature (if any) |
| `genomic_mutations` | optional - metadata for mutations; see below |

The column format changes depending on the `type`:

* If we want to annotate a new aminoacid substitution, then
	* type has to be `mutation`
    * `name` can be any string
    * `site` needs to be a HGVS-format protein variant specification, e.g.
     `p.Ala398Tyr`
    * a new column, `genomic_mutations`, can optionally be added to annotate genomic
    mutations metadata corresponding to the protein mutation. For each
    mutation, we can specify one mor more space-separated single-nucleotide
    substitutions in the HGVS format, preposed by hg19 or hg38 depending on the
    reference genome assembly (hg19/hg38). For instance, `hg38,17:g.7673776G>A`.
    To actually use the information from the 'genomic_mutations' column, the
    `genomic_mutations` metadata should also be specified when calling `ManualAnnotation.add_mutations`
    See the examples below.

* If we want to annotate a post-translational modification, then 
    * `type` should be one of `ptm_phosphorylation`,
`ptm_ubiquitination`, `ptm_acetylation`, `ptm_sumoylation`, `ptm_nitrosylation`,
`ptm_methylation`
    * `name` can be any string
    * `site` needs to be a single number, e.g. `34`. This is the residue number
    in the sequence (1-based) on which the PTM is found
    * `function` can be any string
    * `reference` can be any string

* If we want to annotate a novel short linear motif, then 
	* `type` should be `linear_motif`
	* `site` should be a dash-separated residue range, i.e. `10-25` to signify
	from residue 10 to 25 in the aminoacid sequence, 1-based, including extremities
    * `function` can be any string
    * `reference` can be any string

* If we want to annotate a structured region, then
	* `type` should be `structure`
	* `site` should be a dash-separated residue range, i.e. `10-25` to signify
	from residue 10 to 25 in the aminoacid sequence, 1-based, including extremities
    * `function` can be any string
    * `reference` can be any string

For instance, this is a working example of the csv file (named `test.csv`):

```
name;site;type;function;reference;genomic_mutations
asd;p.Met1Ala;mutation;;qwe;hg38,17:g.7673776G>A
qwe;3;ptm_phosphorylation;asd;
zxc;10-25;linear_motif;zzz;
ert;30-40;structure;xxx;
```

If you do not want to include the genomic_mutations data, then structure the csv file as follows.

```
name;site;type;function;reference
asd;p.Met1Ala;mutation;;qwe
qwe;3;ptm_phosphorylation;asd
zxc;10-25;linear_motif;zzz
ert;30-40;structure;xxx
```

using the CSV file works as you would expect:

```py
>>> ma = ManualAnnotation('test.csv')

# adds mutations to the seq object
>>> ma.add_mutations(seq)
# adds PTM annotations to the sequence object
>>> ma.add_position_properties(seq)

# alternatively, if we want to add genomic mutations metadata, we need to
# explicitly add it in the metadata argument instead. The following commented
# line adds mutations with the corresponding metadata to the seq obejct
# >>> ma.add_mutations(seq, metadata=['genomic_mutations'])

# adds structure or linear motif annotation to the sequence object
>>> ma.add_sequence_properties(seq)
```

We can then verify:

```py

>>> seq.positions[0].mutations
[<Mutation M1A from Manual annotations from test.csv>]

>>> seq.positions[2].properties['ptm_phosphorylation']
<PositionProperty Phosphorylation Site from Manual annotations from test.csv>

>>> seq.properties['linear_motif']
[<SequenceProperty Linear motif from Manual annotations from test.csv, positions 10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25>]

>>> seq.properties['structure']
[<SequenceProperty Structure from Manual annotations from test.csv, positions 30,31,32,33,34,35,36,37,38,39,40>]
```

## Generating the final table

Once the Sequence object has been annotated with the desired data and metadata
we can generate a summary table containing all the information. This is done using
the Table module:

```py
>>> from cancermuts.table import Table

>>> tbl = Table()

# generate pandas data frame
>>> df = tbl.to_dataframe(seq)
```

We can then manipulate the dataframe as we see fit, e.g. by saving it
as a csv file:

```py
# save pandas dataframe as CSV
>>> mt = Table()
>>> df = mt.to_dataframe(seq)
>>> df.to_csv("metatable.csv")
```

## Plotting the content of the final table

Finally, we can use Cancermuts to obtain a graphical representation of the
final data frame:

```py
>>> fig, ax = mt.plot_metatable(df, fname='my_table.pdf', section_size=50)
```

the function returns the figure and axes `Matplotlib` objects relative
to the figure that is being created. It also saves the final output as the file
indicated by the `fname` argument in pdf format. The plot is similar to the one
shown and described in our publication. Different aspects of the plot can
be customized specifying different arguments, the major ones being:

* the `section_size` arguments control how many protein sequence position are
described by each section of the plot, meaning it controls how many sections
are generated. The argument is the number of desired positions
(e.g. `section_size=50`)
* several data types can be hidden or shown by changing options to `False`
(do not show) or `True` (show; default for most cases). These are `mutations`,
`elm`, `ptms`, `structure`, `structure_mobidb`, `mutations_revel`
* by default, only ELMs overlapping mutation sites are plotted - this can be
changed by using argument `mutation_elms_only=False`
* option `figsize` accepts a tuple of number (width and height) and allows to
change size and proportion of the output figure (e.g. `figsize=(4,5)`)

