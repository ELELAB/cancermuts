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
from future.utils import iteritems

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import argparse
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
            'Adaptin binding Endosome-Lysosome-Basolateral sorting signals':'Adaptin',
            'APCC-binding Destruction motifs':'APCC',
            'AP-2 beta2 appendage CCV component motifs':'AP2',
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
            'Cyclin N-terminal Domain Docking Motifs':'Cyclin',
            'Cyclin recognition site':'Cyclin',
            'Di-Tryptophan targeting motif to the Delta-COP MHD domain':'MHD',
            'DDB1-Cullin4 binding site':'DDB1-Cullin4',
            'DLC1/2 binding site':'DLC1/2',
            'Endosome-Lysosome-Basolateral sorting signals':'ELB ss',
            'Extracellular side LRP5 and -6 binding motif':'LRP5,6',
            'EVH1  ligands':'EVH1',
            'FFAT motif':'FFAT',
            'FHA phosphopeptide ligands':'FHA',
            'Helical calmodulin binding motifs': 'Calm.',
            'IAP-binding motif (IBM)':'IBM',
            'Immunoreceptor tyrosine-based motif':'ImmY',
            'Integrin RGD-type binding sites':'Integrin',
            'IRF-3 interaction and dimerisation motif':'IRF-3',
            'IRF-3 binding site':'IRF-3',
            'KEAP1 binding degron':'KEAP1',
            'LC3 binding site':'LC3',
            'LYPxL motif':'LYPxL',
            'MAPK docking motifs':'MAPK',
            'MYND domain binding motif.':'MYND',
            'N-degron':'Ndeg',
            'NES Nuclear Export Signal':'NES',
            'NLS classical Nuclear Localization Signals':'NLS',
            'NRD cleavage site':'NRD cl.',
            'Nuclear receptor box':'NRB',
            'PCNA binding PIP Box':'PCNA',
            'PCSK cleavage site':'PCSK cl.',
            'Phosphotyrosine ligands bound by SH2 domains':'pY',
            'PP1-docking motif RVXF':'PP1',
            'PP1-docking motif SILK':'PP1',
            'PP2A holoenzyme B56-docking site':'PP2A',
            'PP2AC binding site':'PP2AC',
            'PP4 EVH1-binding docking motifs':'PP4 EVH1',
            'PTB ligand':'PTB',
            'PDZ ligands':'PDZ',
            'PDZ domain ligands':'PDZ',
            'Pex14 ligand motif':'Pex14',
            'Raptor interacting motif':'Raptor',
            'RIR motif':'RIR',
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
            'WAVE regulatory complex (WRC) binding site motif':'WRC',
            'WW domain ligands':'WW',
            'WDR5 WD40 repeat (blade 5,6)-binding ligand':'WDR5',
            'WxxL LIR motif':'WxxL',
            'Y-based sorting signal':'',
            'di Arginine retention/retrieving signal':'diArg',
            'eIF4E binding motif':'eIF4E',
            'xLIR LIR motif':'xLIR',
            'GTPase-binding domain (GBD) ligand':'GTPase'}

    ptm_colors = defaultdict(lambda: 'black',
                {   'ptm_acetylation'     : 'grey',
                    'ptm_methylation'     : 'darkgreen',
                    'ptm_ogalnac'         : 'orange',
                    'ptm_oglcnac'          : 'darkorange',
                    'ptm_phosphorylation' : 'red',
                    'ptm_ubiquitination'  : 'blue',
                    'ptm_sumoylation'     : 'lightblue',
                    'ptm_nitrosylation'   : 'cyan',
                    'ptm_cleavage'        : 'magenta'
                        })

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
        ptm_codes = {}

        for k,v in iteritems(position_properties_classes):
            headers[k] = v.header
            if 'ptm' in v.category:
                ptms[k] = v.header
                ptm_codes[k] = v.code

        headers['ptm_sources'] = 'ptm_sources'

        for k,v in iteritems(sequence_properties_classes):
            headers[k] = v.header

        for k,v in iteritems(metadata_classes):
            headers[k] = v.header

        headers['position']    = SequencePosition.header
        headers['wt']          = 'ref_aa'
        headers['mutated']     = 'alt_aa'
        headers['mut_sources'] = 'sources'

        self.headers = headers
        self.ptms = ptms
        self.ptm_codes = ptm_codes

        if not set(self.ptms.keys()).issubset(set(self.ptm_colors.keys())):
            self.log.warning("not enough colors specified for ptms to plot")



    def to_dataframe(self, sequence, mutation_metadata=["cancer_study", "cancer_type", "genomic_coordinates", "genomic_mutations", "revel_score", "cancer_site", "cancer_histology",'gnomad_exome_allele_frequency', 'gnomad_genome_allele_frequency', 
                                                        'gnomad_popmax_exome_allele_frequency', 'gnomad_popmax_genome_allele_frequency', 'clinvar_variant_id', 'clinvar_germline_classification', 'clinvar_germline_condition', 'clinvar_germline_review_status', 
                                                        'clinvar_oncogenicity_classification', 'clinvar_oncogenicity_condition', 'clinvar_oncogenicity_review_status', 'clinvar_clinical_impact_classification', 'clinvar_clinical_impact_condition', 'clinvar_clinical_impact_review_status'],
                        position_properties=['ptm_phosphorylation','ptm_methylation','ptm_ubiquitination','ptm_cleavage', 'ptm_nitrosylation','ptm_acetylation', 'ptm_sumoylation', 'ptm_ogalnac', 'ptm_oglcnac', 'mobidb_disorder_propensity'],
                        sequence_properties=['linear_motif', 'structure']):

        rows = []

        header =  [ self.headers['position'] ]
        header += [ self.headers[p] for p in position_properties ]
        header += [self.headers['ptm_sources']]
        sequence_properties_cols_start = len(header)

        header += [ self.headers[p] for p in sequence_properties ]
        sequence_properties_col = list(range(sequence_properties_cols_start, len(header)))
        header += [self.headers['wt'], self.headers['mutated'], self.headers['mut_sources']]
        for md in mutation_metadata:
            header.append(self.headers[md])

        positions_mutlist = []

        for gi,p in enumerate(sequence.positions):
            base_row = [p.sequence_position]
            ptm_sources_list = []
            for r in position_properties:
                if r in p.properties:
                    val = p.properties[r].get_value_str()
                else:
                    val = None
                base_row.append(val)
            ptm_sources_list = []
            for ptm_key in self.ptms.keys():   
                if ptm_key in p.properties:
                    ptm_sources_list.extend([s.name for s in p.properties[ptm_key].sources])
            
            ptm_sources_str = ",".join(sorted(set(ptm_sources_list))) if ptm_sources_list else None
            base_row.append(ptm_sources_str)
            base_row.extend([None]*len(sequence_properties_col))
            base_row.append(p.wt_residue_type)

            mut_strings = [str(m) for m in p.mutations]
            mut_strings_order = sorted(list(range(len(mut_strings))), key=mut_strings.__getitem__)

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
                this_row.append(",".join([s.name for s in p.mutations[m].sources]))
                for md in mutation_metadata:
                    md_values = []
                    try:
                        self.log.info(f"mutation {p.mutations[m]}")
                        for single_md in p.mutations[m].metadata[md]:
                            if single_md is None:
                                continue
                            this_value = single_md.get_value_str()
                            if this_value is None or this_value != this_value: # hacky check for np.nan
                                continue
                            md_values.append(this_value)
                            self.log.info(f"appending {single_md.get_value_str()}")
                        self.log.info(f"values to be joined {md}: {sorted(list(set(md_values)))}")
                        if "clinvar" not in str(md):
                            md_str = ", ".join(sorted(set(map(str, md_values))))
                        else:
                            md_str = ", ".join(map(str, md_values))
                    except KeyError:
                        md_str = None
                    this_row.append(md_str)
                positions_mutlist.append(gi)
                rows.append(this_row)

        df = pd.DataFrame(rows, columns=header)

        this_col = df.pop(self.headers['wt'])
        df.insert(1, self.headers['wt'], this_col)

        this_col = df.pop(self.headers['mutated'])
        df.insert(2, self.headers['mutated'], this_col)

        positions = df[ self.headers['position'] ]

        property_cols = dict()

        properties_series = dict()

        delenda = []
        for pidx, property_name in enumerate(sequence_properties):
            if property_name not in sequence.properties:
                #property_cols[property_name] = [ None ] * len(positions)
                delenda.append(property_name)
                self.log.warning("property %s not found in data; it won't be annotated" % property_name)
                continue
            elif sequence.properties[property_name][0].category == 'linear_motif':
                property_cols[property_name] = dict( [ (i, list()) for i in list(set(positions)) ] )
            elif sequence.properties[property_name][0].category == 'structure':
                property_cols[property_name] = dict( [ (i, list()) for i in list(set(positions)) ] )

        for d in delenda:
            sequence_properties.remove(d)

        for pidx, property_name in enumerate(sequence_properties):
            props = sequence.properties[property_name]
            for p in props:
                if p.category == 'linear_motif':
                    for pos in p.positions:
                        property_cols[property_name][pos.sequence_position].append(p.get_value_str())
                if p.category == 'structure':
                    for pos in p.positions:
                        property_cols[property_name][pos.sequence_position].append(p.get_value_str())

        for pidx, property_name in enumerate(sequence_properties):
            properties_series[property_name] = []

            if p.category == 'linear_motif':
                for pos in positions:
                    properties_series[property_name].append("|".join(property_cols[property_name][pos]))
                df[self.headers[property_name]] = properties_series[property_name]

            if p.category == 'structure':
                for pos in positions:
                    properties_series[property_name].append("|".join(property_cols[property_name][pos]))
                df[self.headers[property_name]] = properties_series[property_name]

        return df

    def _splice_metatable(self, df, section_size=100):
        # position ranges are left and right-inclusive
        # actual positions are left-inclusive only

        positions = sorted(list(set(df[self.headers['position']])))

        start = positions[0]

        nsegments = len(positions) // section_size

        rest = len(positions) % section_size

        dfs = []
        position_ranges = []

        for i in range(nsegments):
            position_ranges.append((start +  i    * section_size,
                                    start + (i+1) * section_size))
            dfs.append(df[df[self.headers['position']].between(*position_ranges[-1], inclusive='left')])

        if rest > 0:
            position_ranges.append((start + (i+1) * section_size, start + (i+1) * section_size + rest))
            dfs.append(df[df[self.headers['position']].between(*position_ranges[-1], inclusive='left')])
        return dfs, position_ranges

    def _y_ladder(self, miny, maxy, nelm):
        return cycle(np.linspace(miny, maxy, nelm))

    def _plot_mutations(self, ax, df_i, revel, revel_not_annotated, revel_cutoff, y_ladder):

        df_m = df_i[ df_i[self.headers['mutated']].notnull() ]

        if len(df_m) == 0:
            return

        if revel is True:
            def process_revel(x):
                if x == '':
                    return np.nan
                if pd.isna(x) or type(x) is float:
                    return x
                vals = x.split(',')
                return max(map(float, vals))

            pd.options.mode.chained_assignment = None

            numeric_revel = df_m[self.headers['revel_score']].apply(process_revel)
            df_m[self.headers['revel_score']] = numeric_revel
            df_m = df_m.fillna(revel_not_annotated)

            pd.options.mode.chained_assignment = 'warn'

            ax.set_ylabel("REVEL score")
            ax.set_ylim((0.0, 1.0))

        else:
            df_m[self.headers['revel_score']] = revel_not_annotated

        for col in list(set(df_m['stem_colors'])):

            df_m_c = df_m[df_m['stem_colors'] == col]

            scores = []
            for s, row in df_m_c.iterrows():
                try:
                    this_score = float(row[self.headers['revel_score']])
                except ValueError:
                    this_score = np.max(np.array(row[self.headers['revel_score']].split(", "), dtype=float))
                    self.log.warning("The revel score for {0}{1}{2} had more than one value ({3}). The largest will be used".format(
                        row[self.headers['wt']],
                        row[self.headers['position']],
                        row[self.headers['mutated']],
                        row[self.headers['revel_score']]))
                scores.append(this_score)

            (markerline, stemlines, baseline) = ax.stem(df_m_c[self.headers['position']], scores)
            plt.setp(baseline, visible=False)
            plt.setp(stemlines, 'color', col)
            plt.setp(markerline, 'color', col)

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
            if self.ptms[t] not in df_i.keys():
                self.log.warning("PTM type %s not found in DataFrame" % self.ptms[t])
                continue
            df_t = df_i[ df_i[self.ptms[t]].fillna('') == self.ptm_codes[t] ]

            if len(df_t) == 0:
                self.log.info("No information found for %s" % self.ptms[t])
                continue
            for r in df_t[self.headers['position']]:
                ax.text(r, self.y_ptm, self.ptm_codes[t], color=self.ptm_colors[t], horizontalalignment='center') #bbox=dict(facecolor='red', alpha=0.5),
                ax.axvline(x=r, color=self.ptm_colors[t], lw=0.5)

    def _plot_elms(self, ax, df, df_i, mutation_elms_only=True, color='lightblue', color_manual='orange', y_ladder=(-0.2, -0.6, 5)):

        all_elms = []

        df_e = df [ df[ self.headers['linear_motif'] ].notnull() ][ self.headers['linear_motif']]
        df_e = filter(lambda x: x != '', df_e)

        df_mut_pos = set( df[ df[self.headers['mutated']].notnull() ][self.headers['position']] )
        df_i_pos = sorted(list(set(df_i[self.headers['position']])))

        df_i_range = (df_i_pos[0], df_i_pos[-1])

        for e in df_e:
            all_elms.extend(e.split('|'))

        all_elms = set(all_elms)

        ladder = self._y_ladder(*y_ladder)

        all_elms = [ elm.split(", ") for elm in all_elms ]

        for i, elm in enumerate(all_elms):

            if len(elm) < 3:
                raise TypeError(f"format of ELM entry {elm} was unexpected")
            if len(elm) > 3: # this happens if the description contains a comma
                description = ", ".join(elm[0:-2])
                all_elms[i] = [ description, elm[-2], elm[-1] ]
                elm = all_elms[i]

            all_elms[i][1] = [ int(e) for e in elm[1].split('-') ]
            elm_string = all_elms[i][0]
            last_parenthesis_index = elm_string.rfind("(")
            elm_name = elm_string[:last_parenthesis_index].strip()
            elm_code = elm_string[last_parenthesis_index+1:-1]

            elm_first = np.max((df_i_range[0], elm[1][0]))
            elm_last  = np.min((df_i_range[1], elm[1][1]))

            all_elms[i].append( elm_name )
            all_elms[i].append( elm_code )
            all_elms[i].append( elm_first + (elm_last - elm_first) / 2.0 )

        all_elms = sorted(all_elms, key=lambda x: x[3])

        for elm in all_elms:
            _1, pos, source, name, _2, x = elm

            if mutation_elms_only:
                elm_range = set(range(pos[0], pos[1]+1))
                if source != 'Manual annotations':
                    if df_mut_pos.isdisjoint(elm_range):
                        continue

            if pos[0] < df_i_range[0]:
                pos[0] = df_i_range[0]
            if pos[1] > df_i_range[-1]:
                pos[1] = df_i_range[-1]

            # XXX: add colors predicted/manual

            if source == 'Manual annotations':
                this_color = color_manual
            else:
                this_color = color

            ax.add_patch(patches.Rectangle((pos[0],0), pos[1]-pos[0], 1.0, alpha=0.8, color=this_color))

            if x >= df_i_range[0] and x <= df_i_range[-1]:
                if name in self.labels:
                    label = self.labels[name]
                else:
                    label = name
                    self.log.warning("missing ELM label: %s" % name)
                ax.text(x, next(ladder), label, color='black', ha='center', va='center') #fontdict=dict(weight='bold')) #bbox=dict(facecolor='red', alpha=0.5),

    def _plot_structures(self, ax, df, df_i):
        all_structs = []

        df_e = df [ df[ self.headers['structure'] ].notnull() ][ self.headers['structure']]
        df_e = filter(lambda x: x != '', df_e)

        df_i_pos = sorted(list(set(df_i[self.headers['position']])))
        df_i_range = (df_i_pos[0], df_i_pos[-1])

        for e in df_e:
            all_structs.extend(e.split('|'))

        all_structs = set(all_structs)

        all_structs = [ struct.split(", ") for struct in all_structs ]

        for i, struct in enumerate(all_structs):
            all_structs[i][1] = [ int(e) for e in struct[1].split('-') ]
            all_structs[i].append( struct[1][0] + (struct[1][1] - struct[1][0]) / 2.0 )

        all_structs = sorted(all_structs, key=lambda x: x[3])

        for struct in all_structs:
            name, pos, source, x = struct
            if name != 'structured':
                continue
            if pos[0] < df_i_range[0]:
                pos[0] = df_i_range[0]
            if pos[1] > df_i_range[-1]:
                pos[1] = df_i_range[-1]
            ax.add_patch(patches.Rectangle((pos[0],0), pos[1]-pos[0], 1.0, alpha=0.3, color="black", hatch='...', fill=False))

    def _plot_structures_mobidb(self, ax, df, df_i):
        all_structs = []

        df_e = df [ df[ self.headers['mobidb_disorder_propensity'] ].notnull() ][ self.headers['mobidb_disorder_propensity']]
        df_e = filter(lambda x: x != 'S', df_e)

        df_i_pos = sorted(list(set(df_i[self.headers['position']])))
        df_i_range = (df_i_pos[0], df_i_pos[-1])

        from itertools import groupby
        from operator import itemgetter

        ranges =[]

        for k,g in groupby(enumerate(df_i_pos),lambda x:x[0]-x[1]):
            group = (map(itemgetter(1),g))
            group = list(map(int,group))
            ranges.append((group[0],group[-1]))

        all_structs = [(r[0], r[-1]) for r in ranges]

        for pos in all_structs:
            if pos[0] < df_i_range[0]:
                pos[0] = df_i_range[0]
            if pos[1] > df_i_range[-1]:
                pos[1] = df_i_range[-1]
            ax.add_patch(patches.Rectangle((pos[0],0), pos[1]-pos[0], 1.0, alpha=0.3, color="black", hatch='...', fill=False))

    def plot_metatable( self,
                        df,
                        fname=None,
                        section_size=50,
                        figsize=(8.27,6),
                        mutations=True,
                        elm=True,
                        ptms=True,
                        structure=True,
                        structure_mobidb=False,
                        ptm_types=None,
                        mutations_revel=True,
                        revel_not_annotated=0.0,
                        filter_elms=True,
                        y_ladder=(0.65, 0.95, 4),
                        elm_y_ladder=(-0.2, -0.6, 5),
                        revel_cutoff=0.4,
                        stem_colors=None,
                        rcParams={'font.size':8.0, 'font.sans-serif':['Arial']},
                        mutation_elms_only=True):


        if rcParams:
            for k,v in iteritems(rcParams):
                matplotlib.rcParams[k] = v

        df = df.copy()

        stem_colors_present = not stem_colors is None

        if stem_colors_present:
            if "stem_colors" in df.columns:
                self.log.warning("stem_colors column will be overridden")
            if len(stem_colors) != df.shape[0]:
                self.log.error("one stem color per dataframe line needs to be provided")
                raise TypeError

            df['stem_colors'] = stem_colors

        elif not "stem_colors" in df.columns:
            stem_colors = [plt.rcParams['axes.prop_cycle'].by_key()['color'][0]] * df.shape[0]
            self.log.debug("default stem color will be used")

            df['stem_colors'] = stem_colors

        else:
            self.log.debug("stem_colors column will be used for stem colors")

        dfs, position_ranges = self._splice_metatable(df, section_size)

        fig, axes = plt.subplots(nrows=len(dfs), sharex=False, figsize=figsize)
        #                    left=0.02, top=1.02, bottom=0.02, right=1.02)

        if self.headers['mutated'] not in df.keys():
            self.log.warning("No mutations column found - mutations won't be annotated")
            mutations = False

        if self.headers['revel_score'] not in df.keys():
            self.log.warning("No revel score found - value won't be annotated")
            mutations_revel = False

        if self.headers['linear_motif'] not in df.keys():
            self.log.warning("No ELMs information found - won't be annotated")
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
                self._plot_elms(ax, df, df_i, mutation_elms_only=mutation_elms_only, y_ladder=elm_y_ladder)

            if structure:
                self._plot_structures(ax, df, df_i)

            if structure_mobidb:
                self._plot_structures_mobidb(ax, df, df_i)

        plt.tight_layout()
        plt.subplots_adjust(hspace=0.6)


        if fname:
            fig.savefig(fname)

        return fig, axes

    def parse_metatable(self, fname):
        return pd.read_csv(fname, delimiter=';')
