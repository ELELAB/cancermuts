# import the UniProt data source class

from cancermuts.datasources import (
    UniProt, cBioPortal, PhosphoSite, dbPTM, GlyGen, COSMIC,
    MobiDB, MyVariant, RevelDatabase, ManualAnnotation, ClinVar, NetPhos
)
from cancermuts.exceptions import *
from cancermuts.core import ProteinVariant
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

# the seq.sequence_numbering attribute is an ordered list of the protein positions:
print(seq.sequence_numbering[0:3])

# confirm non-canonical status
print("Is the sequence canonical?", seq.is_canonical)
print(
    "Transcript accession for filtering (if any):",
    seq.aliases.get('transcript_accession') or seq.aliases.get('ensembl_transcript_id')
)

# Use COSMIC to retrieve non-canonical alternative isoform annotations
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

cosmic_variants = []
seen = set()
for position in seq.sequence_numbering:
    for variant in seq.variants_at_position(position):
        if variant.hgvs not in seen and any(source.name == "COSMIC" for source in variant.sources):
            cosmic_variants.append((position, variant))
            seen.add(variant.hgvs)
if cosmic_variants:
    first_position, first_mutation = cosmic_variants[0]
    print("First COSMIC mutation:", first_mutation, "at position", first_position)
    print("Sources:", first_mutation.sources)
    print("Metadata:", first_mutation.metadata)
else:
    print(
        "No COSMIC mutations found for transcript:",
        seq.aliases.get('transcript_accession') or seq.aliases.get('ensembl_transcript_id') or "canonical"
    )

# cBioPortal does not support non-canonical mutations
cbioportal = cBioPortal()

try:
    cbioportal.add_mutations(seq)
except UnexpectedIsoformError:
    print("cBioPortal mutations will not be added, as a non-canonical isoform has been provided")

# Use ClinVar to retrieve non-canonical alternative isoform annotations
clinvar = ClinVar()
clinvar.add_mutations(seq, metadata=[
    'clinvar_germline_classification',
    'clinvar_germline_condition',
    'clinvar_germline_review_status',
    'genomic_mutations',
    'clinvar_variant_id',
    'genomic_coordinates',
    'clinvar_oncogenicity_condition',
    'clinvar_oncogenicity_classification',
    'clinvar_oncogenicity_review_status',
    'clinvar_clinical_impact_condition',
    'clinvar_clinical_impact_review_status',
    'clinvar_clinical_impact_classification'
])

clinvar_variants = []
seen = set()
for position in seq.sequence_numbering:
    for variant in seq.variants_at_position(position):
        if variant.hgvs not in seen and variant.metadata.get("clinvar_variant_id"):
            clinvar_variants.append((position, variant))
            seen.add(variant.hgvs)
if clinvar_variants:
    first_position, first_mutation = clinvar_variants[0]
    print("First ClinVar mutation:", first_mutation, "at position", first_position)
    last_position, last_mutation = clinvar_variants[-1]
    if last_position == first_position and last_mutation == first_mutation:
        print("(Only one ClinVar mutation found)")
    else:
        print("Last ClinVar mutation:", last_mutation, "at position", last_position)
else:
    print("No ClinVar mutations found for RefSeq:", seq.aliases["refseq"])

# annotate with REVEL using local database
rl = RevelDatabase("/data/databases/REVEL/revel_with_transcript_ids")
rl.add_metadata(seq)

# print annotated mutation and REVEL score
mut = seq.variants_at_position(820)[0]
print("Mutation:", mut)
print("REVEL score:", mut.metadata.get('revel_score', []))

# PhosphoSite does not support non-canonical isoforms
ps = PhosphoSite('/data/databases/phosphosite/')

try:
    ps.add_sequence_properties(seq)
except UnexpectedIsoformError:
    print("PhosphoSite annotations will not be added, as a non-canonical isoform has been provided")

# dbPTM does not support non-canonical isoforms
db = dbPTM('/data/databases/dbPTM/')

try:
    db.add_sequence_properties(seq)
except UnexpectedIsoformError:
    print("dbPTM annotations will not be added, as a non-canonical isoform has been provided")

# GlyGen does not support non-canonical isoforms
gg = GlyGen('/data/databases/GlyGen/', database_file='human_proteoform_glycosylation_sites_uniprotkb_filtered.csv')

try:
    gg.add_sequence_properties(seq)
except UnexpectedIsoformError:
    print("GlyGen annotations will not be added, as a non-canonical isoform has been provided")

# NetPhos supports non-canonical isoforms through isoform-specific local files
np = NetPhos('/data/databases/netphos_human_proteome/netphos_human_isoforms/raw/')
np.add_sequence_properties(seq)

# MobiDB does not support non-canonical isoforms
mdb = MobiDB()

try:
    mdb.add_sequence_properties(seq)
except UnexpectedIsoformError:
    print("MobiDB annotations will not be added, as a non-canonical isoform has been provided")

# MyVariant does not support non-canonical isoforms
mv = MyVariant()

try:
    mv.add_metadata(seq)
except UnexpectedIsoformError:
    print("REVEL scores from MyVariant annotations will not be added, as a non-canonical isoform has been provided")

# Save MetaTable
tbl = Table()

df = tbl.to_dataframe(seq)
df.to_csv("metatable_non_canonical.csv")
