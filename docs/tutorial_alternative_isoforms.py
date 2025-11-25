# import the UniProt data source class
from cancermuts.datasources import UniProt, cBioPortal, PhosphoSite, COSMIC, MobiDB, MyVariant, RevelDatabase, ManualAnnotation, ClinVar
from cancermuts.exceptions import *
from cancermuts.core import Mutation
from cancermuts.metadata import GenomicMutation
from cancermuts.table import Table

# create the UniProt object
up = UniProt()

# retrieve a specific UniProt isoform
seq = up.get_sequence("AMBRA1", isoform="Q9C0C7-2")

# provide refseq of selected isoform for ClinVar parsing:
seq.aliases["refseq"] = "NP_001287660"

# this prints the downloaded isoform sequence
print(seq.sequence)

# the seq.positions attribute is an ordered list of the protein positions:
print(seq.positions[0:5])

# confirm non-canonical status
print("Is the sequence canonical?", seq.is_canonical)
print("Transcript accession for filtering (if any):",
      seq.aliases.get('transcript_accession') or seq.aliases.get('ensembl_transcript_id'))

# Use COSMIC to retrieve non-canonical alternative isofrom annotations

cosmic = COSMIC(
    targeted_database_file='/data/databases/cosmic-v102/Cosmic_CompleteTargetedScreensMutant_v102_GRCh38.tsv',
    screen_mutant_database_file='/data/databases/cosmic-v102/Cosmic_GenomeScreensMutant_v102_GRCh38.tsv',
    classification_database_file='/data/databases/cosmic-v102/Cosmic_Classification_v102_GRCh38.tsv',
    database_encoding='latin1',
    lazy_load_db=True
)

cosmic.add_mutations(
    seq,
    genome_assembly_version='GRCh38',
    metadata=['genomic_coordinates', 'genomic_mutations', 'cancer_site', 'cancer_histology']
)

positions_with_mut = [p for p in seq.positions if p.mutations]
if positions_with_mut:
    first_pos = positions_with_mut[0]
    first_mut = first_pos.mutations[0]
    print("First COSMIC mutation:", first_mut, "at position", first_pos.sequence_position)
    print("Sources:", first_mut.sources)
    print("Metadata:", first_mut.metadata)
else:
    print("No COSMIC mutations found for transcript:",
          seq.aliases.get('transcript_accession') or seq.aliases.get('ensembl_transcript_id') or "canonical")

# cBioPortal does not support non-canonical mutations
cbioportal = cBioPortal()

try:
    cbioportal.add_mutations(seq)
except UnexpectedIsoformError:
    print("cBioPortal mutations will not be added, as a non-canonical isoform has been provided")

# Use ClinVar to retrieve non-canonical alternative isoform annotations
clinvar = ClinVar()
clinvar.add_mutations(seq, metadata=[
        'clinvar_germline_classification', 'clinvar_germline_condition', 'clinvar_germline_review_status', 'genomic_mutations',
        'clinvar_variant_id', 'genomic_coordinates', 'clinvar_oncogenicity_condition', 'clinvar_oncogenicity_classification',
        'clinvar_oncogenicity_review_status', 'clinvar_clinical_impact_condition', 'clinvar_clinical_impact_review_status', 
        'clinvar_clinical_impact_classification'])

is_cv = lambda m: m.metadata.get("clinvar_variant_id") is not None

positions_with_cv = [
    p for p in seq.positions
    if p.mutations and any(is_cv(m) for m in p.mutations)
]

if positions_with_cv:
    first_pos = positions_with_cv[0]
    first_mut = [m for m in first_pos.mutations if is_cv(m)][0]
    print("First ClinVar mutation:", first_mut, "at position", first_pos.sequence_position)

    last_pos = positions_with_cv[-1]
    last_mutations = [m for m in last_pos.mutations if is_cv(m)]
    last_mut = last_mutations[-1]

    if last_pos == first_pos and last_mut == first_mut:
        print("(Only one ClinVar mutation found)")
    else:
        print("Last ClinVar mutation:", last_mut, "at position", last_pos.sequence_position)
else:
    print("No ClinVar mutations found for refseq:", seq.aliases['refseq'])

#print(seq.positions[819].mutations[0])
#print(seq.positions[819].mutations[0].sources)
#print(seq.positions[819].mutations[0].metadata)

# annotate with REVEL using local database
rl = RevelDatabase("/data/databases/REVEL/revel_with_transcript_ids")
rl.add_metadata(seq)

# print annotated mutation and REVEL score
mut = seq.positions[819].mutations[0]
print("Mutation:", mut)
print("REVEL score:", mut.metadata.get('revel_score', []))

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

# Save MetaTable
tbl = Table()

df = tbl.to_dataframe(seq)
df.to_csv("metatable_non_canonical.csv")
