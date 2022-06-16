from cancermuts.datasources import *
from cancermuts.core import *
from cancermuts.log import *
from cancermuts.table import *

"""
u=UniProt()
s= u.get_sequence('AMRA1_HUMAN')
print s.positions
print len(s.positions)
"""

start_logging()

melanoma_studies = ['mel_tsam_liang_2017', 
                                        'skcm_broad_dfarber', 
                                        'skcm_vanderbilt_mskcc_2015', 
                                        'skcm_broad', 
                                        'skcm_tcga', 
                                        'skcm_yale', 
                                        'skcm_ucla_2016', #?? 
                                        'desm_broad_2015', 
                                        'uvm_tcga'  ] 

melanoma_studies += [  'skcm_tcga_pan_can_atlas_2018',
                       'uvm_tcga_pan_can_atlas_2018' ]

melanoma_studies += ['skcm_mskcc_2014', 'skcm_dfci_2015', 'mel_ucla_2016', 'um_qimr_2016', 'skcm_broad_brafresist_2012']

melanoma_types = ['malignant_melanoma', 'malignant_melanoma_of_soft_parts-clear_cell_sarcoma']

gene_id = "AMBRA1"
gene_id_cosmic = "AMBRA1_ENST00000458649"
gene_id_uniprot = "AMRA1_HUMAN"
up=UniProt()
seq = up.get_sequence(gene_id)
seq.aliases["cosmic"] = gene_id_cosmic
seq.aliases["uniprot"] = gene_id_uniprot

cb=cBioPortal(cancer_studies=melanoma_studies)
cb.add_mutations(seq, metadata=['cancer_type','cancer_study','genomic_coordinates','genomic_mutations'], mainisoform=1)

cosmic_v91=['/data/databases/cosmic-v91/CosmicMutantExport.tsv']
cosmic=COSMIC(database_files=cosmic_v91)
cosmic.add_mutations(seq, cancer_types=melanoma_types, use_alias="cosmic", metadata=['cancer_type', 'genomic_coordinates','genomic_mutations'])

mv=MyVariant()
mv.add_metadata(seq)

gnomad = gnomAD(version='2.1')
gnomad.add_metadata(seq)

ps = PhosphoSite()
ps.add_position_properties(seq)

excludenda = [
"CLV_MEL_PAP_1",
"CLV_PCSK_KEX2_1",
"LIG_BIR_III_4",
"LIG_LIR_Apic_2",
"LIG_LIR_Nem_3",
"LIG_WD40_WDR5_VDV_2",
"CLV_MEL_PAP_1",
"CLV_NRD_NRD_1",
"CLV_PCSK_KEX2_1",
"CLV_PCSK_PC1ET2_1",
"CLV_PCSK_PC1ET2_1",
"CLV_PCSK_SKI1_1",
"CLV_PCSK_PC7_1",
"DEG_SPOP_SBC_1",
"DOC_PP4_FxxP_1",
"DOC_PP4_MxPP_1",
"DOC_USP7_MATH_1",
"DOC_USP7_MATH_2",
"DOC_USP7_UBL2_3",
"LIG_BRCT_BRCA1_1",
"LIG_FHA_1",
"LIG_FHA_2",
"LIG_MYND_1",
"LIG_PCNA_PIPBox_1",
"LIG_PCNA_yPIPBox_3",
"LIG_RGD",
"LIG_SH2_STAP1",
"LIG_SUMO_SIM_anti_2",
"LIG_SUMO_SIM_par_1",
"LIG_TRFH_1",
"LIG_UBA3_1",
"LIG_WD40_WDR5_VDV_2",
"TRG_Pf-PMV_PEXEL_1",
"MOD_.*"]

excludenda_str = '(' + '|'.join(excludenda) + ')'

elm = ELMPredictions()
elm.add_sequence_properties(seq, exclude_elm_classes=excludenda_str, use_alias="uniprot")

gnomad = gnomAD(version='2.1')
gnomad.add_metadata(seq)

ma = ManualAnnotation('manual_annotation.csv', delimiter=';', header=0, index_col=False)
ma.add_position_properties(seq)
ma.add_sequence_properties(seq)
mt = Table()
df = mt.to_dataframe(seq)
df.to_csv("metatable.csv")
