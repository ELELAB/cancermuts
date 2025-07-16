import sys
sys.path.insert(0, "/data/user/elenikia/cancermuts")

# import the UniProt data source class
from cancermuts.datasources import UniProt 


# create the corresponding uniprot object
up = UniProt()

# get the sequence for the protein
seq = up.get_sequence('MAP1LC3B')

# alternatively, we can specifically ask for a Uniprot ID
seq = up.get_sequence('MAP1LC3B', upid='MLP3B_HUMAN')

seq.aliases["refseq"] = "NP_073729"

# this prints the downloaded protein sequence
print(seq.sequence)

# the seq.positions attribute is an ordered list of the protein positions:
seq.positions[0:3]

# import data sources classes
from cancermuts.datasources import cBioPortal, COSMIC, ClinVar

# add mutations from cBioPortal

cb = cBioPortal()

cb.add_mutations(seq, metadata=['cancer_type', 'cancer_study', 'genomic_mutations'])

# let us check out some of the mutations
print(seq.positions[38].mutations)
print(seq.positions[64].mutations)
print(seq.positions[122].mutations)

print(seq.positions[64].mutations[0].sources)

print(seq.positions[64].mutations[0].mutated_residue_type)

print(seq.positions[38].mutations[0].metadata)

# add mutations from COSMIC
cosmic = COSMIC(database_files=['/data/databases/cosmic-v96/CosmicMutantExport.tsv'],
                database_encoding=['latin1'])

cosmic.add_mutations(seq, 
					 metadata=['genomic_coordinates', 'genomic_mutations', 
					 			'cancer_site', 'cancer_histology'])

# let's check them out
print(seq.positions[64].mutations[0])

print(seq.positions[64].mutations[0].sources)

print(seq.positions[64].mutations[0].metadata)

# add mutations from ClinVar
clinvar = ClinVar()
clinvar.add_mutations(seq, metadata=[
    'clinvar_classification',
    'clinvar_condition',
    'clinvar_review_status',
    'clinvar_genomic_annotation',
    'clinvar_method',
    'clinvar_variant_id',
    'clinvar_variant_name'
])

# Check ClinVar metadata (variant ID, condition, classification, etc.)
print(seq.positions[15].mutations[0].metadata.get('clinvar_variant_id'))
print(seq.positions[15].mutations[0].metadata.get('condition'))
print(seq.positions[15].mutations[0].metadata.get('classification'))
                                        
# add annotations from MyVariant (REVEL)
from cancermuts.datasources import MyVariant

mv = MyVariant()
mv.add_metadata(seq)

print(seq.positions[64].mutations[0].metadata['revel_score'])

# add annotations from gnomAD
from cancermuts.datasources import gnomAD

gnomad = gnomAD(version='2.1')
gnomad.add_metadata(seq, md_type=['gnomad_exome_allele_frequency',
	                              'gnomad_genome_allele_frequency'])

print(seq.positions[64].mutations[0].metadata['gnomad_exome_allele_frequency'])

print(seq.positions[64].mutations[0].metadata['gnomad_genome_allele_frequency'])

from cancermuts.datasources import PhosphoSite, MobiDB

# add annotations from PhosphoSite
ps = PhosphoSite('/data/databases/phosphosite/')
ps.add_position_properties(seq)

print(seq.positions[4].properties)

print(seq.positions[28].properties)

# add annotations from MobiDB
mdb = MobiDB()
mdb.add_position_properties(seq)

print(seq.positions[0].properties['mobidb_disorder_propensity'])

print(seq.positions[10].properties['mobidb_disorder_propensity'])

# add annotations from ELM
from cancermuts.datasources import ELMPredictions

elm = ELMPredictions()
elm.add_sequence_properties(seq,
			    exclude_elm_classes="MOD_.")

print(seq.properties)
print(seq.properties['linear_motif'][0].type)

# save table
from cancermuts.table import Table

tbl = Table()

df = tbl.to_dataframe(seq)

df.to_csv("metatable_pancancer.csv")

tbl.plot_metatable(df, fname='my_table_pancancer.pdf', section_size=50)
