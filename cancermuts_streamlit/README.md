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

### Downloading database

No database currently available, an example database is available for testing.

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
streamlit run main_page.py
```

a web browser window pointing to the web app should open
