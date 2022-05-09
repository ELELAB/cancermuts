## Installation

Cancermuts is easily installed as a regular Python package on any supported
platform, as follows:

1. clone the cancermuts clonebase from GitHub to a local folder of your choice.

```
cd my_cancermuts
git clone https://github.com/ELELAB/cancermuts.git
```

2. in order to install it, we recommend creating a virtual environment for it
so that it is isolated from the rest of the system. For instance, using
virtualenv package, if available:

```
virtualenv -p python3.7 ./venv
```

another popular choice would be using a conda environment:

```
conda create cancermuts python=3.7
```

3. install the package. This should also install all the requirements:

```
python setup.py install
```

4. you're now ready to use Cancermuts! Please head to the tutorial section 

{% hint style='info' %}
You can use cancermuts by means of a Python script or interactively, in a
[Jupyter notebook](https://jupyter.org) or [iPython](https://ipython.org).
{% endhint %}

## Offline resources

As it will become more clear later, Cancermuts expect some resources to be
available on the local computer to be used. Currently, this is only true
for the COSMIC Mutant Export, which contains a list of curated cancer mutations.

In order to obtain it
1. head to the [COSMIC website](https://cancer.sanger.ac.uk/cosmic), 
register and log in

2. head to the Data menu at the top and select the Downloads option

3. head to the "COSMIC Mutation Data" section and downlaod the corresponding
Whole file (CosmicMutantExport.tsv.gz)

4. decompress the file, e.g. using `gzip`:
```
gunzip -d CosmicMutantExport.tsv.gz
```

{% hint style='warning' %}
The downloaded file will take a significant amount of disk space, at least 15GB.
{% endhint %}

{% hint style='tip' %}
Make note of the location of this file as it will come handy in the Tutorial
{% endhint %}
