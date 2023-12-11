# Cancermuts - Streamlit application
# Copyright (C) 2023 Matteo Tiberti, Danish Cancer Society
#               2023 Elena Papaleo, Danish Cancer Society 
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import streamlit as st
from streamlit_utils import add_affiliation_logo

st.set_page_config(layout="wide",
    page_title="Acknowledgements and data usage",
    page_icon="🙏")

add_affiliation_logo()

st.header('Acknowledgements and data usage')

st.subheader("Citing us")

st.write("""Cancermuts is developed by the Cancer Structural Biology lab at the
Danish Cancer Institute, Copenhagen, Denmark.""")

st.write("""If you use data from this database in your research, please cite:

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

st.subheader("Data and software availability")

st.markdown('''
We also would like to acknowledge our data sources. The datasets we have used 
are available under a variety of licensing terms, which are listed below for 
reuse.

Unless stated otherwise in the following text, Cancermuts data is released under
the [Creative Commons Attribution 4.0 International (CC BY 4.0) license](https://creativecommons.org/licenses/by/4.0/).

Bulk download of Cancermuts data tables is available at [our OSF repository](https://osf.io/)

The Cancermuts software and source code for this website is available at our
[Cancermuts GitHub repository](https://github.com/ELELAB/Cancermuts)
and released under open source license.''')

st.subheader('Mutations')

st.markdown('''Cancermuts uses the following datasets for the mutation data:

  - [**COSMIC**](https://cancer.sanger.ac.uk/cosmic): The Catalogue Of Somatic
    Mutations In Cancer. We use version 96 in our current releases. The data in
    Cancermuts is released according to the COSMIC Non-Commercial Terms and
    Conditions, and released in Cancermuts for non-commercial use only.
  - [**cBioPortal**](https://www.cbioportal.org/) for cancer genomics. Data from
    cBioPortal present in Cancermuts is released under the [Open Data Commons Open
    Database License (ODbL) v1.0](https://opendatacommons.org/licenses/odbl/1-0/)

By default, following the MAVISp protocol, we also include mutations found
in ClinVar:

  - [**ClinVar**](https://www.ncbi.nlm.nih.gov/clinvar/), the NIH database of
    clinically relevant genetic variants. It is released under public domain.''')

st.image("static/NCBI_powered.png")

st.subheader('Mutations metadata')

st.markdown('''Cancermuts uses the following datasets for the mutation metadata:

   - [**REVEL**](https://sites.google.com/site/revelgenomics/) scores were
   downloaded by [myvariants.info](https://myvariant.info/), whose license
   terms are available on their website. REVEL scores are available as
   integrated by dbSNP which is released under public domain.
   - [**gnomAD**](https://gnomad.broadinstitute.org/), database of human genetic
    variation. It is released under the public domain.''')

st.subheader('Protein metadata and predictions')

st.markdown('''Cancermuts uses the following datasets or software for the prediction or annotation of protein-level features, or structure:

   - [**UniProt**](https://www.uniprot.org). It is released under the Creative
  Commons Attribution 4.0 International (CC BY 4.0) License
  - the [**Protein Data Bank**](https://www.rcsb.org). It is released under
  public domain.
  - [**MobiDB**](https://mobidb.org/) for the prediction of protein disorder.
  It is released under the Creative Commons Attribution 4.0 International (CC
  BY 4.0 DEED) License.
  - [**ELM**](http://elm.eu.org), the Eukaryotic Linear Motif resource for functional sites in proteins.
  - [**PhosphoSitePlus**](https://www.phosphosite.org), a database of
  post-translational modifications. Data from PhosphoSitePlus in Cancermuts, (i.e.
  the PTMs column in our database files) is released according to the
  [PhosphoSitePlus' terms and conditions](https://www.phosphosite.org/staticDownloads), and it is not available for commercial use.''')
