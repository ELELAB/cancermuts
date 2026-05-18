# import the UniProt data source class
from cancermuts.datasources import UniProt

# create the corresponding uniprot object
up = UniProt()

# get the sequence for the protein
seq = up.get_sequence('MAP1LC3B')

# alternatively, we can specifically ask for a Uniprot ID
seq = up.get_sequence('MAP1LC3B', upid='MLP3B_HUMAN')

# OPTIONAL: Use a specific UniProt isoform instead of canonical
# seq = up.get_sequence("AMBRA1", isoform='Q9C0C7-2')

# Specify the desired RefSeq isoform for ClinVar parsing:
seq.aliases["refseq"] = "NP_073729"

# this prints the downloaded protein sequence
print(seq.sequence)

# the seq.sequence_numbering attribute is an ordered list of the protein positions:
print(seq.sequence_numbering[0:3])

# import data sources classes
from cancermuts.datasources import cBioPortal, COSMIC, ClinVar

# add mutations from cBioPortal
cb = cBioPortal(cancer_studies=[
    'coadread_dfci_2016',
    'coadread_genentech',
    'coadread_tcga_pan_can_atlas_2018'
])

cb.add_mutations(seq, metadata=['cancer_type', 'cancer_study', 'genomic_mutations'])

# let us check out some of the mutations
print(seq.variants_at_position(39))
print(seq.variants_at_position(65))
print(seq.variants_at_position(123))

print(seq.variants_at_position(65)[0].sources)
print(seq.variants_at_position(65)[0].alt)
print(seq.variants_at_position(39)[0].metadata)

# add mutations from COSMIC
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
    cancer_sites=['large_intestine'],
    cancer_site_subtype_1=['colon'],
    cancer_types=['carcinoma'],
    cancer_histology_subtype_1=['adenocarcinoma'],
    metadata=['genomic_coordinates', 'genomic_mutations', 'cancer_site', 'cancer_histology']
)

# let's check them out
print(seq.variants_at_position(65)[0])
print(seq.variants_at_position(65)[0].sources)
print(seq.variants_at_position(65)[0].metadata)

# add mutations from ClinVar
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

# Check ClinVar Variant
print(seq.variants_at_position(15))
print(seq.variants_at_position(15)[0].sources)
print(seq.variants_at_position(15)[0].alt)
print(seq.variants_at_position(15)[0].metadata['clinvar_germline_condition'])
print(seq.variants_at_position(15)[0].metadata['clinvar_germline_classification'])

# add annotations from MyVariant (REVEL)
from cancermuts.datasources import MyVariant

mv = MyVariant()
mv.add_metadata(seq)

print(seq.variants_at_position(65)[0].metadata['revel_score'])

# add annotations from gnomAD
from cancermuts.datasources import gnomAD

gnomad = gnomAD(version='2.1')
gnomad.add_metadata(seq, md_type=[
    'gnomad_exome_allele_frequency',
    'gnomad_genome_allele_frequency',
    'gnomad_popmax_exome_allele_frequency',
    'gnomad_popmax_genome_allele_frequency'
])

print(seq.variants_at_position(65)[0].metadata['gnomad_exome_allele_frequency'])
print(seq.variants_at_position(65)[0].metadata['gnomad_genome_allele_frequency'])

from cancermuts.datasources import PhosphoSite, dbPTM, GlyGen, MobiDB, NetPhos

# add annotations from PhosphoSite
ps = PhosphoSite('/data/databases/phosphosite/')
ps.add_sequence_properties(seq)

print(seq.properties_at_position(5))
print(seq.properties_at_position(29))

# add annotations from dbPTM
db = dbPTM('/data/databases/dbPTM/')
db.add_sequence_properties(seq)

print(seq.properties_at_position(19, 'ptm_glycosylation'))
print(seq.properties_at_position(29, 'ptm_phosphorylation'))

# add annotations from GlyGen
gg = GlyGen(
    '/data/databases/GlyGen/',
    database_file='human_proteoform_glycosylation_sites_uniprotkb_filtered.csv'
)
gg.add_sequence_properties(seq)
print(seq.properties_at_position(46, 'ptm_glycosylation'))

# add annotations from NetPhos
np = NetPhos('/data/databases/netphos_human_proteome/netphos_human_isoforms/raw/')
np.add_sequence_properties(seq)
print(seq.properties_at_position(3, 'ptm_phosphorylation'))

# add annotations from MobiDB
mdb = MobiDB()
mdb.add_sequence_properties(seq)

print(seq.properties_at_position(1,'mobidb_disorder_propensity'))
print(seq.properties_at_position(11,'mobidb_disorder_propensity'))

# add annotations from ELM
from cancermuts.datasources import ELMPredictions

elm = ELMPredictions()
elm.add_sequence_properties(
    seq,
    exclude_elm_classes="MOD_."
)

print(seq.properties)
print(seq.properties['linear_motif'][0].type)

# save table
from cancermuts.table import Table

tbl = Table()

df = tbl.to_dataframe(seq)

df.to_csv("metatable.csv")

tbl.plot_metatable(df, fname='my_table.pdf', section_size=50)
