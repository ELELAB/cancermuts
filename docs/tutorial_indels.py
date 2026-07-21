# import the UniProt data source class
from cancermuts.datasources import UniProt

# create the corresponding uniprot object
up = UniProt()

# get the sequence for the protein
seq = up.get_sequence('RAD51C')

# alternatively, we can specifically ask for a Uniprot ID
seq = up.get_sequence('RAD51C', upid='RA51C_HUMAN')

# Specify the desired RefSeq isoform for ClinVar parsing:
seq.aliases["refseq"] = "NP_478123"

# this prints the downloaded protein sequence
print(seq.sequence)

# the seq.sequence_numbering attribute is an ordered list of the protein positions:
print(seq.sequence_numbering[0:3])

# import data sources classes
from cancermuts.datasources import cBioPortal, COSMIC, ClinVar

# define the variant types we want to retrieve
indel_types = ("deletion", "insertion", "delins")

# add mutations from cBioPortal
cb = cBioPortal()

cb.add_mutations(seq,
                 metadata=['cancer_type', 'cancer_study', 'genomic_mutations'],
                 variant_types=indel_types)

# let us check out some of the mutations
print(seq.variants_at_position(176)) 
print(seq.variants_at_position(184)) 
print(seq.variants_at_position(176)[0].sources) 
print(seq.variants_at_position(176)[0].alt) 
print(seq.variants_at_position(176)[0].metadata)

# add mutations from COSMIC
cosmic = COSMIC(targeted_database_file='/data/databases/cosmic-v102/Cosmic_CompleteTargetedScreensMutant_v102_GRCh38.tsv',
                screen_mutant_database_file='/data/databases/cosmic-v102/Cosmic_GenomeScreensMutant_v102_GRCh38.tsv',
                classification_database_file='/data/databases/cosmic-v102/Cosmic_Classification_v102_GRCh38.tsv',
                database_encoding='latin1',
                lazy_load_db=True)

cosmic.add_mutations(seq,
                     genome_assembly_version='GRCh38',
                     metadata=['genomic_coordinates', 'genomic_mutations', 'cancer_site', 'cancer_histology'],
                     variant_types=indel_types)

# let's check them out
print(seq.variants_at_position(38))  
print(seq.variants_at_position(176)) 
print(seq.variants_at_position(38)[0]) 
print(seq.variants_at_position(38)[0].sources) 
print(seq.variants_at_position(38)[0].metadata['genomic_mutations'])

# add mutations from ClinVar
clinvar = ClinVar()
clinvar.add_mutations(seq,
                      metadata=['clinvar_germline_classification',
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
                                'clinvar_clinical_impact_classification'],
                     variant_types=indel_types)

# Check ClinVar variants
print(seq.variants_at_position(28))
print(seq.variants_at_position(338))
print(seq.variants_at_position(6)[0].sources)
print(seq.variants_at_position(6)[0].alt)
print(seq.variants_at_position(6)[0].metadata['clinvar_germline_condition'])

# add annotations from gnomAD
from cancermuts.datasources import gnomAD

gnomad = gnomAD(version='2.1')
gnomad.add_metadata(seq, md_type=['gnomad_exome_allele_frequency',
                                  'gnomad_genome_allele_frequency',
                                  'gnomad_popmax_exome_allele_frequency',
                                  'gnomad_popmax_genome_allele_frequency'])

print(seq.variants_at_position(338)[0].metadata['gnomad_exome_allele_frequency'])

# save table
from cancermuts.table import Table

tbl = Table()

df = tbl.to_dataframe(seq, hgvsp=True)

df.to_csv("metatable_indels.csv")

tbl.plot_metatable(df, fname='my_table_indels.pdf', section_size=50)

