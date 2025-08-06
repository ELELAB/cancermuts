# import the UniProt data source class
from cancermuts.datasources import UniProt, cBioPortal, PhosphoSite, COSMIC. MobiDB
from cancermuts.exceptions import *

# create the UniProt object
up = UniProt()

# retrieve a specific UniProt isoform
seq = up.get_sequence("AMBRA1", isoform="Q9C0C7-2")

# this prints the downloaded isoform sequence
print(seq.sequence)

# the seq.positions attribute is an ordered list of the protein positions:
print(seq.positions[0:5])

# confirm non-canonical status
print("Is the sequence canonical?", seq.is_canonical)

# cBioPortal does not support non-canonical isoforms
cbioportal = cBioPortal()

try:
    cbioportal.add_mutations(seq)
except UnexpectedIsoformError:
    print("cBioPortal mutations will not be added, as a non-canonical isoform has been provided")

cosmic = COSMIC(targeted_database_file='/data/databases/cosmic-v102/Cosmic_CompleteTargetedScreensMutant_v102_GRCh38.tsv',
                screen_mutant_database_file='/data/databases/cosmic-v102/Cosmic_GenomeScreensMutant_v102_GRCh38.tsv',
                classification_database_file='/data/databases/cosmic-v102/Cosmic_Classification_v102_GRCh38.tsv',
                transcript_database_file='/data/databases/cosmic-v102/Cosmic_Transcripts_v102_GRCh38.tsv',
                database_encoding='latin1', lazy_load_db=True)

# COSMIC does not support non-canonical isoforms
try:
    cosmic.add_mutations(seq)
except UnexpectedIsoformError:
    print("COSMIC mutations will not be added, as a non-canonical isoform has been provided")

# PhosphoSite does not support non-canonical isoforms
ps = PhosphoSite('/data/databases/phosphosite/')

try:
    ps.add_position_properties(seq)
except UnexpectedIsoformError:
    print("PhosphoSite annotations will not be added, as a non-canonical isoform has been provided")

# MobiDB does not suport non-canonical isoforms
mdb = MobiDB()

try:
    mdb.add_position_properties(seq)
except UnexpectedIsoformError:
    print("MobiDB annotations will not be added, as a non-canonical isoform has been provided")
