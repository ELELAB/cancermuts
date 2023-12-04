# main_page.py for cancermuts
# (c) 2023 Alberte Heering Estad <ahestad@outlook.com>
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

import streamlit as st
from streamlit_utils import add_affiliation_logo

st.set_page_config(layout="wide",
    page_title="Welcome",
    page_icon="👋")

add_affiliation_logo()

st.header("Cancermuts")

st.subheader("Introduction")

st.write("""Cancermuts is a Python package for the retrieval and annotation of
cancer mutations and their context. This database contains results obtained using
the Cancermuts package in the context of [MAVISp](https://www.biorxiv.org/content/10.1101/2022.10.22.513328v4),
our framework for the prediction of cancer mutations, following our standard
protocol for MAVISp.

Cancermuts is developed at the Danish Cancer Institute, and is available as a
Python package at our [GitHub repository](https://www.github.com/ELELAB/cancermuts).

While this website contains results obtained using Cancermuts, the package itself
is not required to access the data. We encourage interested users to use the Cancermuts
package directly to retrieve data and annotations more tailored to their specific
case of study.

If you use these data in your research, please cite:

  - our [Cancermuts paper](https://www.nature.com/articles/s41419-022-05318-2)

  > Tiberti, M., Di Leo, L., Vistesen, M.V. et al. The Cancermuts software package
  > for the prioritization of missense cancer variants: a case study of AMBRA1 in
  > melanoma. Cell Death Dis 13, 872 (2022). 
  > https://doi.org/10.1038/s41419-022-05318-2

  - our [MAVISp preprint](https://www.biorxiv.org/content/10.1101/2022.10.22.513328v4.abstract):

  > Matteo Arnaudi, Ludovica Beltrame, Kristine Degn, Mattia Utichi, Simone Scrima,
  > Pablo Sánchez-Izquierdo Besora, Karolina Krzesińska, Alberte Heering Estad, 
  > Francesca Maselli, Terézia Dorčaková, Jordan Safer, Katrine Meldgård, 
  > Philipp Becker, Valentina Sora, Alberto Pettenella, Julie Bruun Brockhoff, 
  > Amalie Drud Nielsen, Jérémy Vinhas, Peter Wad Sackett, Claudia Cava, Anna Rohlin, 
  > Mef Nilbert, Sumaiya Iqbal, Matteo Lambrughi, Matteo Tiberti, Elena Papaleo
  > bioRxiv 2022.10.22.513328; doi: https://doi.org/10.1101/2022.10.22.513328
""")

st.subheader("Accessing the data")

st.write("""Please see the Database section on the left for the available data.
The data is available in the form of a CSV file, which can be downloaded either
from this website or from our [OSF repository](https://osf.io/jc32x/)""")

st.subheader("Data usage")

st.write("""Please see the Acknowledgments section on the left for data usage policy and 
further citation information.""")
