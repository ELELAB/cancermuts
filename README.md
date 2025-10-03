Cancer Structural Biology, Danish Cancer Society Research Center, 2100, Copenhagen, Denmark

Repository associated to the publication:

The Cancermuts software package for the prioritization of missense cancer variants: a case study of AMBRA1 in melanoma. 
Matteo Tiberti\*, Luca Di Leo, Mette Vixø Vistesen, Rikke Kuhre, Francesco Cecconi, Daniela De Zio\*, Elena Papaleo.
Submitted to bioRxiv. https://doi.org/10.1101/2022.05.23.493014

contacts for this repository: elenap-at-cancer.dk, tiberti-at-cancer.dk

## cancermuts

### Description

Cancermuts is a package for the automtic retrieval and annotation of cancer
mutations. At the moment it supports:

* Retrieval of protein sequence from UniProt
* Annotation of the sequence with data from
    * cBioPortal - cancer mutations
    * COSMIC - cancer mutations (after SNP filtering)
    * ClinVar - disease mutations
    * ELM - linear motifs
    * REVEL - pathogenicity prediction
    * Phosphosite - post-translational modifications
    * gnomAD - allele frequency
    * MobiDB - protein disorder
    * Manually curated annotations from csv files

Please see the [Cancermuts documentation on GitLab](https://matteo-tiberti.gitbook.io/cancermuts/) for further details

### Cancermuts web server

We have collected a number of Cancermuts entry in connection with the [MAVISp project](https://github.com/ELELAB/mavisp)
and others. We are releasing this dataset on a webserver, hosted by the Technical University of Denmark,
Department of Health Technologies. The webserver includes the possibility to download raw data for each
entry as well as some graphical visualization of the results, that can be tweaked in real time.

[The webserver is available at this link](https://services.healthtech.dtu.dk/services/Cancermuts-1.0/)

If you use the data from the webserver, please cite both the Cancermuts and the MAVISp paper, as 
described in the Acknowledgements section of the website.

