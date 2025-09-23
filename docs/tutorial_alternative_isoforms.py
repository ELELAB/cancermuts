# import the UniProt data source class
from cancermuts.datasources import UniProt, cBioPortal, PhosphoSite, COSMIC, MobiDB, MyVariant, RevelDatabase
from cancermuts.exceptions import *
from cancermuts.core import Mutation
from cancermuts.metadata import GenomicMutation

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

# define mutation manually: Val820Met on isoform 2
wt = "V"
mut = "M"
position_number = 820

# find position object
site_seq_idx = seq.seq2index(position_number)
position = seq.positions[site_seq_idx]

# create Mutation object
mutation = Mutation(position, mut, sources=["Manual"])

# add genomic mutation metadata manually
gm = GenomicMutation(source="Manual", genome_build="hg19", definition="11:g.46456582C>T")
mutation.metadata['genomic_mutations'] = [gm]

# add mutation to position
position.add_mutation(mutation)

# annotate with REVEL using local database
rl = RevelDatabase("/data/databases/REVEL/revel_with_transcript_ids")
rl.add_metadata(seq)

# print annotated mutation and REVEL score
mut = seq.positions[site_seq_idx].mutations[0]
print("Mutation:", mut)
print("REVEL score:", mut.metadata.get('revel_score', []))

# cBioPortal does not support non-canonical mutations
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

# MyVariant does not suport non-canonical isoforms
mv = MyVariant()

try:
    mv.add_metadata(seq)
except UnexpectedIsoformError:
    print("REVEL scores from MyVariant annotations will not be added, as a non-canonical isoform has been provided")



