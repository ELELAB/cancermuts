# import the UniProt data source class
from cancermuts.datasources import UniProt, cBioPortal, PhosphoSite, MobiDB
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

# cBioPortal does not support non-canonical mutations
cbioportal = cBioPortal()

try:
    cbioportal.add_mutations(seq)
except UnexpectedIsoformError:
    print("cBioPortal mutations will not be added, as a non-canonical isoform has been provided")

# PhosphoSite does not support non-canonical mutations
ps = PhosphoSite('/data/databases/phosphosite/')

try:
    ps.add_position_properties(seq)
except UnexpectedIsoformError:
    print("PhosphoSite annotations will not be added, as a non-canonical isoform has been provided")

# MobiDB does not suport non-canonical positions
mdb = MobiDB()

try:
    mdb.add_position_properties(seq)
except UnexpectedIsoformError:
    print("MobiDB annotations will not be added, as a non-canonical isoform has been provided")
