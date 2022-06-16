# coding: utf-8
import pandas
import pandas as pd
t = pd.read_csv('metatable.csv')
t
t = pd.read_csv('metatable.csv', index=0)
t = pd.read_csv('metatable.csv', index_col=0)
t
t.columns()
t.columns
t['Mutated residue']
t[not t['Mutated residue'].is_nan()]
t[not t['Mutated residue'] is str]
t[not pd.isnull(t['Mutated residue'])]
pd.isnull(t['Mutated residue'])
t[not pd.isnull(t['Mutated residue'])]
not t
get_ipython().run_line_magic('pinfo', ' pd.isnull')
t[~pd.isnull(t['Mutated residue'])]
t['Revel score'] >= 0.4
t[t['Revel score'] >= 0.4]
t[t['Revel score'] < 0.4]
t[~pd.isnull(t['Mutated residue'])]
muts = t[~pd.isnull(t['Mutated residue'])]
muts
muts['Revel score']
muts[pd.isnull(muts.revel_score)]
muts[pd.isnull(muts['Revel score'])]
no_revel_muts = muts[pd.isnull(muts['Revel score'])]
no_revel_muts['Genomic mutation']
no_revel_muts['Genomic coordinates']
no_revel_muts
t[t['Revel score'] < 0.4]
t[t['Revel score'] >= 0.4]
28/70
28.0/70.9
28.0/70.0
t
t.columns
structured = t[~pd.isnull(t['Structure'])]
unstructured = t[pd.isnull(t['Structure'])]
structured
unstructured
structured[~pd.isnull(t['Mutated residue'])]
unstructured[~pd.isnull(t['Mutated residue'])]
structured[~pd.isnull(t['Mutated residue'])]
unstructured[~pd.isnull(t['Mutated residue'])]
unstructured[(~pd.isnull(t['Mutated residue']))]
structured[~structured.isnull(t['Mutated residue'])]
structured[~structured.isnull(structured['Mutated residue'])]
structured[~pd.isnull(structured['Mutated residue'])]
unstructured[~pd.isnull(unstructured['Mutated residue'])]
structured[(~pd.isnull(structured['Mutated residue']))]
structured[(~pd.isnull(structured['Mutated residue'])) & (structured['Revel score'] >=0.4)]
unstructured[(~pd.isnull(unstructured['Mutated residue'])) & (unstructured['Revel score'] >=0.4)]
structures
structured
structured
muts
muts['WT residue']
muts['WT residue']
muts.columns
muts['WT residue', 'Mutated residue']
muts[['WT residue', 'Mutated residue']]
muts[['WT residue', 'Position', 'Mutated residue']]
print(muts[['WT residue', 'Position', 'Mutated residue']])
with pd.option_context('display.max_rows', None, 'display.max_columns', None): print(muts[['WT residue', 'Position', 'Mutated residue']])
muts['WT residue'].value_counts()
muts['Mutated residue'].value_counts()
get_ipython().run_line_magic('history', '')
damaging = muts[muts['Revel score'] >= 0.4]
non_damaging = muts[muts['Revel score'] < 0.4]
damaging['WT residue'].countt_values()
damaging['WT residue'].count_values()
damaging['WT residue'].value_counts()
damaging['Mutated residue'].value_counts()
muts
muts[['WT residue', 'Mutated residue']]
muts['WT residue']
muts['WT residue'] + muts['Mutated residue']
muts['WT residue'] + muts['Mutated residue'].count_values()
muts['WT residue'] + muts['Mutated residue'].values_count()
muts['WT residue'] + muts['Mutated residue'].value_counts()
(muts['WT residue'] + muts['Mutated residue']).value_counts()
structured['WT residue']
structured['WT residue']
muts
structured
get_ipython().run_line_magic('save', 'explore.py 1-88')
