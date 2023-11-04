
# Cancermuts web app

## Introduction

Cancermuts is a software that retrieves information on known cancer-associated mutations from a variety of data sources. This application 
is made to make exploration and visualization of Cancermuts data convenient and accessible to provide a solid preliminary insight into 
the effects of certain mutations.

## Requirements

This application was developed with Python 3.8.0 with the following packages installed:

- streamlit
- pandas
- matplotlib
- upsetplot
- numpy
- cancermuts

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
pip install streamlit pandas matplotlib upsetplot numpy
```
4. For the cancermuts package you'll have to clone the Cancermuts GitHub repository to a local folder:
```
git clone https://github.com/ELELAB/cancermuts.git
```
5. Then install the package to your virtual environment:
```
pip install ./cancermuts
```
### Downloading database

No database currently available, an example database is available for testing and exploration.

### Running the app

To run the app locally you can run the following command from within the cancermuts_streamlit folder of the cloned repository:
```
streamlit run main_page.py
```

