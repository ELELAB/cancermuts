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
on our GitHub repository so we can help you.

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
# (these are usually equivalent)
seq = up.get_sequence('MAP1LC3B', upid='MLP3B_HUMAN')

# this prints the downloaded protein sequence
print(seq.sequence)

# the seq.positions attribute is an ordered list of the protein positions:
seq.positions[0:3]
```


### Collecting cancer mutations

Now that we have the Sequence ready, we can start adding annotations to it.
We will first download cancer mutations from cBioPortal and COSMIC.

For cBioPortal, we first create the respective source object:

```py 
# import data sources classes
from cancermuts.datasources import cBioPortal, COSMIC

cb = cBioPortal()
```

In this case, cBioPortal will gather mutations from all cancer studies available
in cBioPortal.

It is possible to specify a set of cancer studies Cancermuts should use instead,
if we are interested in e.g. a certain cancer type or some other specific studies, 
using the cBioPortal study identifier:

```py
cb = cBioPortal(cancer_studies=['brca_broad', 'brca_sander'])
```

Finally, we use this object to gather the mutations from the selected cancer
studies and add them to our Sequence object we created earlier:

```py
cb.
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

