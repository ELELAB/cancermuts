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

## Tutorial steps

### The Sequence object

Cancermuts works 
The first operation to be performed when starting to work with Cancermuts is creating
a Sequence object. This object represents the protein sequence we want to consider for
the analysis and links together all the information we are going to collect. In order to
create it we use a data source class for UniProt from which we can downlaod protein sequences.
This is only source for sequences available at the moment, more can be added in the future.

To download the main UniProt sequence, you will need the corresponding HGCN symbol only.
The UniProt object automatically identifies the most likely UniProt ID corresponding to the
gene name. Alternatively, it is possible to provide the UniProt ID manually as an option.

This is done as follows:


```py
# import the UniProt data source class
from cancermuts.datasources import UniProt

# create the corresponding uniprot object
up = UniProt()

# get the sequence for the protein
seq = up.get_sequence('MAP1LC3B')

# alternatively, we can specifically ask for a Uniprot ID
seq = up.get_sequence('MAP1LC3B', upid='MLP3B_HUMAN')

# this prints the downloaded protein sequence
print(seq.sequence)

# the seq.positions attribute is an ordered list of the protein positions:
seq.positions[0:3]
```

### Collecting cancer mutations

Now that we have the Sequence ready, we can start adding annotations to it.
We will first download cancer mutations from cBioPortal and COSMIC.

We first import the required datasource classes for them:

```py
# import data sources classes
from cancermuts.datasources import cBioPortal, COSMIC
```

#### cBioPortal

For cBioPortal, we first create the respective source object:

```py
cb = cBioPortal()
```

In this case, cBioPortal will gather mutations from all cancer studies available
in cBioPortal.

It is possible to specify a set of cancer studies Cancermuts should use instead,
if we are interested in e.g. a certain cancer type or some other specific studies, 
using the cBioPortal study identifier:

```py
cb = cBioPortal(cancer_studies=['coadread_dfci_2016', 
	                            'coadread_genentech',
	                            'coadread_tcga_pan_can_atlas_2018'])
```

{% hint style='info' %}
Finding the correct set of identifier for the studies of interest can
be tricky, as it isn't obvious from the cBioPortal website which study
IDs are connected to specific studies. When creating a cBioPortal object
without the `cancer_studies` argument, two attributes of the object can
be helpful to this regard. `cb.cancer_types` is a list of cancer types
as classified by cBioPortal and `cb.cancer_studies` is the full list of
studies available on cBioPortal. These list are regular pandas dataframes
and can be manipulated as such, including saved to csv files or similar.

Alternatively, one can visit the (https://www.cbioportal.org/datasets)[cBioPortal datasets page]
which contains a list of all the available studies. While this list doesn't
contain a column with study IDs, the web links that link the studies in the
"Name" column do contain the IDs. For instance, the study named 
"Breast Invasive Carcinoma (Broad, Nature 2012)" links the page
https://www.cbioportal.org/study?id=brca_broad - meaning the corresponding
study ID is `brca_broad`.
{% endhint %}

Finally, we use this object to gather the mutations from the selected cancer
studies and add them to our Sequence object we created earlier:

```py
cb.add_mutations(seq)
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
cb.add_mutations(seq, metadata=['cancer_type', 'cancer_study', 'genomic_mutations'])
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
In [27]: print(seq.positions[38].mutations)
[<Mutation K39R from cBioPortal>]

In [28]: print(seq.positions[64].mutations)
[<Mutation K65E from cBioPortal>]

In [29]: print(seq.positions[122].mutations)
[<Mutation L123S from cBioPortal>]
```

Each mutation is recorded in a Mutation object that can be further explored:

```py
In [31]: seq.positions[64].mutations
[<Mutation K65E from cBioPortal>]

In [9]: seq.positions[64].mutations[0].sources
[<cancermuts.datasources.cBioPortal at 0x1510a0455c50>]

In [32]: seq.positions[64].mutations[0].mutated_residue_type
'E'

In [12]: seq.positions[38].mutations[0].metadata
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
cosmic = COSMIC(database_files=['/data/databases/cosmic-v95/CosmicMutantExport.tsv'])
```

here the `database_files` argument is a list of strings, each of them is a
database file to be considered. Usually, the argument for this file would be
the COSMIC database file that was downloaded as detailed in the Install section.

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

Here we restrict the search to those mutations that are involved in colon
adenocarcinoma. It should be noted that COSMIC supports up to four 
cancer histology types and four cancer histology subtypes and Cancermuts
allows to filter by any of this. In particular, only the mutations that 
correspond to *all* the selected criteria are retained. Please see the
(https://cancer.sanger.ac.uk/cosmic/classification)[COSMIC phenotype classification]
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
In [13]: seq.positions[64].mutations[0]
<Mutation K65E from cBioPortal,COSMIC>

In [14]: seq.positions[64].mutations[0].sources
[<cancermuts.datasources.cBioPortal at 0x1510a0455c50>,
 <cancermuts.datasources.COSMIC at 0x15109ff10050>]

In [15]: seq.positions[64].mutations[0].metadata
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
(https://myvariant.info)[MyVariant.info] database. REVEL is associated to
genomic mutations and not to amino-acid mutations, therefore our amino acid
mutations need to be annotated with the genomic mutations metadata in order
for this to be possible. This is performed as detailed above. In order to
download the associated REVEL score, we can use the appropriate data source
class:

```py
from cancermuts.datasources import MyVariant

mv = MyVariant()
mv.add_metadata(seq)
```

we can check that the mutations have been annotated with the REVEL score:

```py
seq.positions[64].mutations[0].metadata['revel_score']
Out[23]: [<DbnsfpRevel, 0.314>, <DbnsfpRevel, 0.314>]
```

In this case we see the same score twice, as the mutation was previously
annotated with two genomic mutations. The two genomic mutations corresponded
to the same mutations annotated in two different genome assemblies, therefore
the two scores we are able to gather have the same value.

#### gnomAD allele frequencies

Similarly, we annotate mutations with their exome or genome allele frequencies
as found in the (https://www.gnomad.org)[gnomAD database]. This works as you 
would expect by now:

```py
from cancermuts.datasources import gnomAD

gnomad = gnomAD(version='2.1')
gnomad.add_metadata(seq, md_type=['gnomad_exome_allele_frequency',
	                              'gnomad_genome_allele_frequency'])
```

here, we specify the `version` argument to specify the version of gnomAD to be 
considered. Please refer to the API documentation for all the available versions.

The `md_type` keyword allows to select which metadata type(s) to annotate, i.e.
choose between exome or allele frequency (or both as in the example).

We calculate the allele frequency as the ratio between the total allele count over
the allele number as found in gnomAD, if the entry for the corresponding variant 
is available.

Finally we can check the downloaded metadata:

```py
In [31]: seq.positions[64].mutations[0].metadata['gnomad_exome_allele_frequency']
[<gnomADExomeAlleleFrequency, 0.000028>,
 <gnomADExomeAlleleFrequency, 0.000028>]

In [32]: seq.positions[64].mutations[0].metadata['gnomad_genome_allele_frequency']
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
ps = PhosphoSite('/data/databases/phosphosite/')
```

by default, Cancermuts expect the file names in the database to be the default
ones. However it's possible to specify the file name for each post-translational
modifications by supplying the `database_files` argument, which should be a
dictionary associating each PTM to a file name:

```py
my_databse_files = {  'acetylation'     : 'my_Acetylation_site_dataset',
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
ps.add_position_properties(seq, 
	                       properties=['phosphorylation', 'ubiquitination'])
```

However, in this case we will consider all of them:

```py
ps.add_position_properties(seq)
```

Once again, we can check the result:

```py
>>> seq.positions[4].properties
{'ptm_ubiquitination': <PositionProperty Ubiquitination Site from PhosphoSite>}

>>> seq.positions[28].properties
{'ptm_phosphorylation': <PositionProperty Phosphorylation Site from PhosphoSite>}
```

#### structured regions with MobiDB

We use the MobiDB website to predict structured or unstructured regions in our
protein of interest. This works similarly as before:

```py
mdb = MobiDB()
mdb.add_position_properties(seq)
```

we can then check the annotation as done previously:

```py
>>> seq.positions[0].properties['mobidb_disorder_propensity']
<StructuralDisorder, C>

>>> seq.positions[10].properties['mobidb_disorder_propensity']
<StructuralDisorder, S>
```

here an annotation of "C" means coil (disoredered), while an annotation of "S"
means "structured".

### Sequence properties

We further annotate properties to the sequence, i.e. properties that cover
multiple residues. We first import the respective classes:

```py
from cancermuts.datasources import ELMPredictions
```

#### Predictions of short linear motifs from ELM

We annotate the sequence using predictions for short linear motifs from the
[Eukaryotic Linear Motif](http://elm.eu.org) database:

```py
elm = ELMPredictions()
elm.add_sequence_properties(seq,
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

The column format changes depending on the `type`:

* If we want to annotate a new aminoacid substitution, then
	* type has to be `mutation`
    * `name` can be any string
    * `site` needs to be a HGVS-format protein variant specification, e.g.
     `p.Ala398Tyr`
    * `function` should be either empty or a HGVS-format single-nucleotide
    substitution in HGVS format, preposed by 19 or 38 depending on the
    reference genome assembly (hg19/38). For instance, `38,17:g.7673776G>A`

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
name;site;type;function;reference
asd;p.Met1Ala;mutation;38,17:g.7673776G>A;qwe
qwe;3;ptm_phosphorylation;asd;qwe
zxc;10-25;linear_motif;zzz;qqq
ert;30-40;structure;xxx;ppp
```

using the CSV file works as you would expect:

```py
ma = ManualAnnotation('test.csv')

# adds mutations to the seq object
ma.add_mutations(seq )

# adds PTM annotations to the sequence object
ma.add_position_properties(seq )

# adds structure or linear motif annotation to the sequence object
ma.add_sequence_properties(seq )
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
from cancermuts.table import Table

tbl = Table()

# generate pandas data frame
df = tbl.to_dataframe(seq)
```

We can then manipulate the dataframe as we see fit, e.g. by saving it
as a csv file:

```py
# save pandas dataframe as CSV
mt = Table()
df = mt.to_dataframe(seq)
df.to_csv("metatable.csv")
```


