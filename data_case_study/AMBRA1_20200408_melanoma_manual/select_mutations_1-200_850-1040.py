# coding: utf-8
import pandas as pd

data = pd.read_csv('metatable.csv', index_col=0)

muts = data[(data['Position'].between(1, 200) | data['Position'].between(850, 1040)) & ~ data['Mutated residue'].isna()]

muts['chain'] = 'A'

selected_muts = muts['WT residue'] + muts['chain'] + muts['Position'].astype(str) + muts['Mutated residue']

selected_muts.to_csv('mutations_1-200_850-1040.txt', sep=';', index=False, header=False)
