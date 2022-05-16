# Installation

## Requirements

Cancermuts is a Python package and requires at minimum Python 3.7 to run. All
modern linux distributions and Mac OS include a distribution of Python that you
can start using right away.

Furthermore, Cancermuts relies on a number of open source packages, such as:

* requests
* bioservices
* myvariant
* pyliftover
* Bio
* bravado
* matplotlib
* pandas
* parse
* urllib3
* future

Please refer to the `setup.py` file in the distribution for an up to date
list of requirements. It should be noted that requirements are installed
automatically during the installation process described below.

The installation steps also assume you have either the virtualenv Python package
 or Conda installed on your system. You can install the virtualenv package
 by running:

```
$ python -m pip install --user virtualenv
```

that should give you access to the `virtualenv` executable:

```
$ virtualenv
Running virtualenv with interpreter /usr/bin/python2
You must provide a DEST_DIR
Usage: virtualenv.py [OPTIONS] DEST_DIR
...
```

if you'd rather use the Conda package manager, you can follow the instructions
available on their [website](https://docs.conda.io/en/latest/miniconda.html)

The following instructions are meant for Mac OS or a Linux-based operating
 system. Cancermuts should also run on Windows, even though
it hasn't been tested on it and we currently don't support it. If you want to 
use Cancermuts on Windows, we recommend doing so through the [Windows Subsystem
for Linux](https://docs.microsoft.com/en-us/windows/wsl/install) or by means
of a container (e.g. [Docker](https://www.docker.com))

## Installation procedure

Cancermuts is easily installed as a regular Python package on any supported
platform, as follows:

1. create a local copy of the Cancermuts GitHub repository, in a local folder
of your choice:

    ```
    mkdir my_cancermuts && cd my_cancermuts
    git clone https://github.com/ELELAB/cancermuts.git
    ```

2. in order to install it, we recommend creating a virtual environment for it
so that it is isolated from the rest of the system. For instance, using
virtualenv package, if available:

    ```
    virtualenv -p python3.7 ./cancermuts_env
    ```

    another popular choice would be using a conda environment:

    ```
    conda create cancermuts_env python=3.7
    ```

3. Activate the Python environment, for instance if you used virtualenv:

    ```
    source cancermuts_env/bin/activate
    ```

    or conda:

    ```
    conda activate ./cancermuts_env
    ```

    either way you should end up with an active Python environment, i.e.
    there should be a `(cancermuts_env)` string before your command prompt:

    ```
    teo@my-computer:my_folder$ . cancermuts_env/bin/activate
    (cancermuts_env) teo@my-computer:my_folder$
    ```

4. install the package. This should also install all the requirements:

    ```
    pip install ./cancermuts
    ```

5. you're now ready to use Cancermuts!

{% hint style='info' %}
You can use cancermuts by means of a Python script or interactively, in a
[Jupyter notebook](https://jupyter.org) or [iPython](https://ipython.org).
{% endhint %}

## Offline resources

As it will become more clear later, Cancermuts expects some resources to be
available on the local computer to be used. Currently, this is only true
for the COSMIC Mutant Export, which contains a list of cancer mutations.

In order to obtain it:

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

