# table.py - cancermuts table plotting and saving
# (c) 2019 Matteo Tiberti <matteo.tiberti@gmail.com>
# This file is part of cancermuts
#
# cancermuts is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Nome-Programma is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Nome-Programma.  If not, see <http://www.gnu.org/licenses/>.


"""
datasources classes --- :mod:`cancermuts.table`
================================================================
Classes to save/read the table and generate a graphical representation

"""

from __future__ import division

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import argparse
import logging as log
from .log import *
from .properties import *
from .metadata import *
from .core import *
import numbers
from itertools import cycle
from collections import defaultdict
import numpy as np

class Table:
    labels={"'WD' motif binding TPR of kinesin light chain":'WD',
            '14-3-3 binding phosphopeptide motif':'14-3-3',
            'APCC activator-binding ABBA motif':'ABBA',
            'Actin-binding motifs':'Actin',
            'Atg8 protein family ligands':'Atg8',
            'BCL-2 binding site':'BCL-2',
            'BRCT phosphopeptide ligands':'BRCT',
            'Binding motif for UBA3 adenylation domain':'UBA3',
            'COP1 E3 ligase binding degron motif':'COP1 E3',
            'Calcineurin (PP2B)-docking motif LxvP':'LxvP',
            'Caspase cleavage motif':'Caspase cl.',
            'Cks1 ligand':'CKS1',
            'Clathrin box':'Clathrin',
            'CtBP ligand motif':'CtBP',
            'Cyclin docking motif':'Cyclin',
            'Cyclin recognition site':'Cyclin',
            'DDB1-Cullin4 binding site':'DDB1-Cullin4',
            'DLC1/2 binding site':'DLC1/2',
            'Endosome-Lysosome-Basolateral sorting signals':'ELB ss',
            'Extracellular side LRP5 and -6 binding motif':'LRP5,6',
            'EVH1  ligands':'EVH1',
            'FHA phosphopeptide ligands':'FHA',
            'IAP-binding motif (IBM)':'IBM',
            'Immunoreceptor tyrosine-based motif':'ImmY',
            'IRF-3 binding site':'IRF-3',
            'KEAP1 binding degron':'KEAP1',
            'LC3 binding site':'LC3',
            'MAPK docking motifs':'MAPK',
            'MYND domain binding motif.':'MYND',
            'N-degron':'Ndeg',
            'NES Nuclear Export Signal':'NES',
            'NLS classical Nuclear Localization Signals':'NLS',
            'NRD cleavage site':'NRD cl.',
            'Nuclear receptor box':'NRB',
            'PCSK cleavage site':'PCSK cl.',
            'PP1-docking motif RVXF':'PP1',
            'PP2A holoenzyme B56-docking site':'PP2A',
            'PP2AC binding site':'PP2AC',
            'PTB ligand':'PTB',
            'PDZ ligands':'PDZ',
            'Pex14 ligand motif':'Pex14',
            'SCF ubiquitin ligase binding Phosphodegrons':'SCF',
            'Separase cleavage motif':'Separase cl.',
            'SH2 ligand':'SH2',
            'SH3 ligand':'SH3',
            'Siah binding Motif':'Siah',
            'SPOP SBC docking motif':'SPOP',
            'SUMO interaction site':'SUMO int',
            'TRAF2 binding site':'TRAF2',
            'TRAF6 binding site':'TRAF6',
            'Tankyrase-binding motif':'Tankyrase',
            'UEV Domain binding PTAP motif':'UEV',
            'USP7 binding motif':'USP7',
            'WW domain ligands':'WW',
            'WxxL LIR motif':'WxxL',
            'Y-based sorting signal':'',
            'di Arginine retention/retrieving signal':'diArg',
            'eIF4E binding motif':'eIF4E',
            'xLIR LIR motif':'xLIR',
            'GTPase-binding domain (GBD) ligand':'GTPase'}

    ptms = {    'Phosphorylation site'  : 'P',
                'Methylation site'      : 'Me',
                'Ubiquitination site'   : 'Ubq',
                'Caspase cleavage site' : 'C',
                'S-Nitrosylation site'  : 'SN',
                'Acetylation site'      : 'Ac',
                'Sumoylation site'      : 'Sumo'    }

    ptm_colors = defaultdict(lambda: 'black',
                {    'Phosphorylation site'  : 'red',
                    'Methylation site'      : 'darkgreen',
                    'Ubiquitination site'   : 'blue',
                    'Caspase cleavage site' : 'purple',
                    'S-Nitrosylation site'  : 'orange',
                    'Acetylation site'      : 'grey',
                    'Sumoylation site'      : 'lightblue'    })

    headers = { "elm"      : "Residue part of a linear motif",
                "position" : 'Position',
                "mutated"  : 'Mutated residue',
                "wt"       : 'WT residue',
                "revel"    : 'Revel score' }

    y_ptm = 1.02

    @logger_init
    def __init__(self, labels=None, ptm_colors=None, y_ptm=1.02):
        if labels:
            self.labels = labels

        if ptm_colors:
            self.ptm_colors = ptm_colors

        self.y_ptm = y_ptm

        headers = {}
        ptms = {}

        for k,v in position_properties_classes.iteritems():
            headers[k] = v.description
            if 'ptm' in v.category:
                ptms[k] = v.description

        for k,v in sequence_properties_classes.iteritems():
            headers[k] = v.description

        for k,v in metadata_classes.iteritems():
            headers[k] = v.description

        headers['position']    = SequencePosition.description
        headers['wt']          = 'Wt residue'
        headers['mutated']     = 'Mutated residue'
        headers['mut_sources'] = 'Mutation sources'

        self.headers = headers

        if not set(self.ptms.keys()).issubset(set(self.ptm_colors.keys())):
            self.log.warning("not enough colors specified for ptms to plot")



    def to_dataframe(self, sequence, mutation_metadata=["cancer_study", "cancer_type", "genomic_coordinates", "genomic_mutations", "revel_score", "cancer_site", "cancer_histology"], 
                        position_properties=['ptm_phosphorylation','ptm_methylation','ptm_ubiquitination','ptm_cleavage', 'ptm_nitrosylation','ptm_acetylation', 'ptm_sumoylation','mobidb_disorder_propensity'],
                        sequence_properties=['linear_motif',]):

        rows = []

        header =  [ self.headers['position'] ]
        header += [ self.headers[p] for p in position_properties ]
        sequence_properties_cols_start = len(header)

        header += [ self.headers[p] for p in sequence_properties ]
        sequence_properties_col = range(sequence_properties_cols_start, len(header))
        header += ['WT residue', 'Mutated residue', 'Sources']
        for md in mutation_metadata:
            header.append(self.headers[md])

        positions_mutlist = []

        for gi,p in enumerate(sequence.positions):
            base_row = [p.sequence_position]
            print p.properties
            for r in position_properties:
                if r in p.properties.keys():
                    print "addng", p.properties.keys(), r
                    val = p.properties[r].get_value_str()
                else:
                    print "notin", p.properties.keys(), r
                    val = None
                base_row.append(val)
            base_row.extend([None]*len(sequence_properties_col))
            base_row.append(p.wt_residue_type)

            mut_strings = [str(m) for m in p.mutations]
            mut_strings_order = sorted(range(len(mut_strings)), key=mut_strings.__getitem__)
            
            if len(mut_strings_order) == 0:
                this_row = list(base_row)
                this_row.append(None)
                this_row.append(None)
                this_row.extend([None]*len(mutation_metadata))
                positions_mutlist.append(gi)
                rows.append(this_row)
                continue

            for m in mut_strings_order:
                this_row = list(base_row)
                this_row.append(p.mutations[m].mutated_residue_type)
                this_row.append(",".join( [s.name for s in p.mutations[m].sources] ))
                for md in mutation_metadata:
                    md_values = []
                    try:
                        for single_md in p.mutations[m].metadata[md]:
                            md_values.append(single_md.get_value_str())

                        md_str = ", ".join(sorted(list(set(md_values))))
                    except:
                        md_str = None
                    this_row.append(md_str)
                #print "appending", gi
                positions_mutlist.append(gi)
                #print positions_mutlist
            rows.append(this_row)

        print rows
        print header

        df = pd.DataFrame(rows, columns=header)

        print df

        positions = df[ self.headers['position'] ]

        property_cols = dict()

        properties_series = dict()

        for pidx, property_name in enumerate(sequence_properties):
            if property_name not in sequence.properties.keys():
                property_cols[property_name] = [ None ] * len(positions)
                sequence_properties.remove(property_name)
                self.log.warning("property %s not found in data; it won't be annotated" % property_name)
                continue
            elif sequence.properties[property_name][0].category == 'linear_motif':
                property_cols[property_name] = dict( [ (i, list()) for i in list(set(positions)) ] )

        for pidx, property_name in enumerate(sequence_properties):
            props = sequence.properties[property_name]
            for p in props:
                if p.category == 'linear_motif':
                    for pos in p.positions:
                        property_cols[property_name][pos.sequence_position].append(p.get_value_str())

        for pidx, property_name in enumerate(sequence_properties):
            if p.category == 'linear_motif':
                properties_series[property_name] = []
                for pos in positions:
                    properties_series[property_name].append("|".join(property_cols[property_name][pos]))
            df[self.headers[property_name]] = properties_series[property_name]

        return df

    def _splice_metatable(self, df, section_size=100):

        positions = sorted(list(set(df[self.headers['position']])))

        nsegments = len(positions) // section_size

        rest = len(positions) % section_size

        dfs = []
        position_ranges = []

        for i in range(nsegments):
            position_ranges.append((i * section_size + 1, (i+1) * section_size + 1))
            dfs.append(df[df[self.headers['position']].between(*position_ranges[-1])])

        if rest > 0:
            position_ranges.append(((i+1) * section_size + 1, (i+1) * section_size + 1 + rest))
            dfs.append(df[df[self.headers['position']].between(*position_ranges[-1])])
        return dfs, position_ranges

    def _y_ladder(self, miny, maxy, nelm):
        return cycle(np.linspace(miny, maxy, nelm))


    def _plot_mutations(self, ax, df_i, revel, revel_not_annotated, revel_cutoff, y_ladder):

        df_m = df_i[ df_i[self.headers['mutated']].notnull() ]

        if revel is True:
            df_m = df_m.fillna(revel_not_annotated)
            ax.set_ylabel("REVEL score")
            ax.set_ylim((0.0, 1.0))

        else:
            df_m['Revel score'] = revel_not_annotated

        (markerline, stemlines, baseline) = ax.stem(df_m[self.headers['position']], df_m[self.headers['revel']])
        plt.setp(baseline, visible=False)

        ladder = self._y_ladder(*y_ladder)

        all_positions = sorted(list(set(df_m[self.headers['position']])))
        for p in all_positions:
            this_muts = df_m[df_m[self.headers['position']] == p]
            this_muts_str = ", ".join(["%s%d%s" % tuple(v) for v in this_muts[[self.headers['wt'], self.headers['position'], self.headers['mutated']]].values])
            ax.text(p, next(ladder), this_muts_str, horizontalalignment='center', verticalalignment='center')

        if revel_cutoff:
            ax.axhline(y=revel_cutoff, color='gray', lw=0.5, ls='--')


    def _plot_ptms(self, ax, df_i, types):
        if types is None:
            types = self.ptms.keys()

        for t in types:
            df_t = df_i[ df_i[t] == self.ptms[t] ]
            for r in df_t[self.headers['position']]:
                ax.text(r, self.y_ptm, self.ptms[t], color=self.ptm_colors[t], horizontalalignment='center') #bbox=dict(facecolor='red', alpha=0.5),
                ax.axvline(x=r, color=self.ptm_colors[t], lw=0.5)

    def _plot_elms(self, ax, df, df_i, mutation_elms_only=True, color='lightblue', y_ladder=(-0.2, -0.6, 5)):
        
        all_elms = []

        df_e = df [ df[ self.headers['elm'] ].notnull() ][ self.headers['elm']]

        df_mut_pos = set( df[ df[self.headers['mutated']].notnull() ][self.headers['position']] )
        df_i_pos = sorted(list(set(df_i[self.headers['position']])))
        df_i_range = (df_i_pos[0], df_i_pos[-1])

        for e in df_e:
            all_elms.extend(e.split('|'))

        all_elms = set(all_elms)

        ladder = self._y_ladder(*y_ladder)

        all_elms =     [ elm.split(", ") for elm in all_elms ]
        for i, elm in enumerate(all_elms):
            all_elms[i][2] = map(int, elm[1].split('-'))
            all_elms[i].append( elm[1][0] + (elm[1][1] - elm[1][0]) / 2.0 )

        all_elms = sorted(all_elms, key=lambda x: x[4])

        for elm in all_elms:
            name, pos, source, x = elm
            if pos[0] < df_i_range[0]:
                pos[0] = df_i_range[0]
            if pos[1] > df_i_range[-1]:
                pos[1] = df_i_range[-1]

            # XXX: add colors predicted/manual
            if mutation_elms_only:
                elm_range = set(range(pos[0], pos[1]+1))
                if df_mut_pos.isdisjoint(elm_range):
                    continue

            ax.add_patch(patches.Rectangle((pos[0],0), pos[1]-pos[0], 1.0, alpha=0.8, color=color))

            if x >= df_i_range[0] and x <= df_i_range[-1]:
                if name in self.labels.keys():
                    label = self.labels[name]
                else:
                    label = name
                    self.log.warning("missing ELM label: %s" % name)
                ax.text(x, next(ladder), label, color='black', ha='center', va='center') #fontdict=dict(weight='bold')) #bbox=dict(facecolor='red', alpha=0.5),

    def plot_metatable( self, 
                        df, 
                        fname=None, 
                        section_size=50, 
                        figsize=(8.27,6), 
                        mutations=True, 
                        ptms=True, 
                        elm=True, 
                        ptm_types=None, 
                        mutations_revel=True, 
                        revel_not_annotated=0.5, 
                        do_elms=True, 
                        filter_elms=True, 
                        y_ladder=(0.65, 0.95, 4), 
                        revel_cutoff=0.4,
                        rcParams={'font.size':8.0, 'font.sans-serif':['Arial']}):

        if rcParams:
            for k,v in rcParams.iteritems():
                matplotlib.rcParams[k] = v

        df = df.copy()

        dfs, position_ranges = self._splice_metatable(df, section_size)

        fig, axes = plt.subplots(nrows=len(dfs), sharex=False, figsize=figsize)
        #                    left=0.02, top=1.02, bottom=0.02, right=1.02)

        if self.headers['mutated'] not in df.keys():
            #log.warning("No mutations column found - mutations won't be annotated")
            mutations = False

        if self.headers['revel'] not in df.keys():
            #log.warning("No revel score found - value won't be annotated")
            mutations_revel = False

        if self.headers["elm"] not in df.keys():
            #log.warning("No ELMs information found - won't be annotated")
            elm = False

        for i, df_i in enumerate(dfs):

            ax = axes[i]
            ax.set_xlim((position_ranges[i][0], position_ranges[i][1]-1))
            ax.set_ylim((0.0, 1.0))
            ax.minorticks_on()
            ax.tick_params(axis='y',which='minor',bottom='off')

            if mutations:
                self._plot_mutations(ax, df_i, mutations_revel, revel_not_annotated, revel_cutoff, y_ladder)

            if ptms:
                self._plot_ptms(ax, df_i, ptm_types)

            if elm:
                self._plot_elms(ax, df, df_i)

        plt.tight_layout()
        plt.subplots_adjust(hspace=0.6)


        if fname:
            fig.savefig(fname)

        return fig, axes

    def parse_metatable(self, fname):
        return pd.read_csv(fname, delimiter=';')
