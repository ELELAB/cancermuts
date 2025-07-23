# import the UniProt data source class
from cancermuts.datasources import UniProt, cBioPortal
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

