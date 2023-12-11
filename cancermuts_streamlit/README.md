# Cancermuts web app

## Introduction

Cancermuts is a software that retrieves information on known cancer-associated
mutations from a variety of data sources. This web app is made to make
exploration and visualization of Cancermuts data convenient and accessible
to provide a solid preliminary insight into the effects of certain mutations.

## Requirements

In principle, it is compatible with all operating systems that support Python.
It requires Python >=3.8 and the following Python packages:

- streamlit 1.28.2
- pandas 2.1.3
- matplotlib 3.8.2
- upsetplot 0.8.0
- numpy 1.26.2
- cancermuts, from commit ID 4795011
- bioservices 1.11.2
- osfclient 0.0.5

It has been last test on Linux (Ubuntu 18.04), and on macOS (13.5.2),
with Python 3.8.0 and the following package versions:

This application was last tested with Python 3.9.6 and the following packages installed:

- streamlit 1.28.2
- pandas 2.1.3
- matplotlib 3.8.2
- upsetplot 0.8.0
- numpy 1.26.2
- cancermuts, commit ID 4795011
- bioservices 1.11.2
- osfclient 0.0.5

## Usage

### Installation

It is recommended to use the Python environment module `virtualenv` to manage the required packages.
`
1. Create virtual environment:
```
virtualenv -p python cancermuts_env
```
2. Activate it:
```
source cancermuts_env/bin/activate
```
3. Install requirements using `pip`:
```
pip install streamlit==1.28.2 pandas==2.1.3 matplotlib==3.8.2 upsetplot==0.8.0 numpy==1.26.2 bioservices==1.11.2
```
4. Install the cancermuts package, by cloning the Cancermuts GitHub repository to a local folder and then
installing the package in your virtual environment:
```
git clone https://github.com/ELELAB/cancermuts.git
cd cancermuts
pip install .
```

### Running the app - full database

In order to run our web server locally with its full content, you will need to
download the full Cancermuts dataset from OSF, as follow, as well as download our
web app from GitHub. If you'd rather test the web app on a small subset, please
follow the instructions in the "Running the app - test dataset" instead.

1. If you haven't already, activate your Python environment (see previous steps)

2. enter your local copy of the Cancermuts repository in your system, and the
`cancermuts_streamlit` subfolder:

```
cd cancermuts/cancermuts_streamlit
```

3. Download the database files from [our OSF repository](https://osf.io/jc32x/).
If you have installed the `osfclient` Python package (see requirements), just run:

```
rm -rf ./database
osf clone && mv jc32x/osfstorage/database/ . && rm -r jc32x
```

alternatively, you can download the whole OSF database from the web interface
and copy its `database` folder in `cancermuts_streamlit`. At the end of the process,
you should have the OSF `database` folder and its contents inside the `cancermuts_streamlit` folder.

4. With your Python environment still active and from inside the `cancermuts_streamlit` repository
directory, run:

```
streamlit run Welcome.py
```

a browser window displaying the MAVISp web app should open.

### Running the app - test database

After installing all requirements, in order to run the web app you will need
to have your virtualenv Python environment still active (see previous instructions).
If this is the case,

1. Enter the `cancermuts_streamlit` folder. We will need to set a systen variable
to make the web app aware of the location of the test database:

```
cd cancermuts_streamlit
export CANCERMUTS_DATABASE=./example_data/example_database
```

2. start the web app:

```
streamlit run Welcome.py
```

a web browser window pointing to the web app should open
