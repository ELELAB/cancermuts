import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('Agg')
from cancermuts.datasources import *
from cancermuts.core import *
from cancermuts.log import *
from cancermuts.table import *
import pandas as pd
"""
u=UniProt()
s= u.get_sequence('AMRA1_HUMAN')
print s.positions
print len(s.positions)
"""
start_logging()

mt = Table()
df = pd.read_csv('metatable.csv')

af_co = 0.0001   # 0.01%

def_color = plt.rcParams['axes.prop_cycle'].by_key()['color'][0]
#stem_colors = [def_color] * df.shape[0]
#df['stem_colors'] = [ 'grey' if c > af_co else def_color for c in df['Exome allele frequency (gnomAD v2.1)'] ]

mt.plot_metatable(df, "plot.pdf", section_size=300, figsize=(8.27,13.0), y_ladder=(0.55,0.95,5), elm_y_ladder=(-0.2, -0.5, 4))


