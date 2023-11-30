# main_page.py for cancermuts
# (c) 2023 Alberte Heering Estad <ahestad@outlook.com>
# This file is part of cancermuts
#
# cancermuts is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cancermuts is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with cancermuts.  If not, see <http://www.gnu.org/licenses/>.

import streamlit as st
import os
import pandas as pd
from streamlit_utils import *
import matplotlib.pyplot as plt
from upsetplot import plot, from_memberships
import numpy as np
import io
from cancermuts.table import Table

def download_plot_button(filename):
    fn = filename+'.png'
    img = io.BytesIO()
    plt.savefig(img, format='png')
    
    btn = st.download_button(
                label="Download plot",
                data=img,
                file_name=fn,
                mime="image/png")
    
def filter_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    # addFilter = st.checkbox("Add filter")
    
    df['Date'] = pd.to_datetime(df['Date'])

    filter_container = st.container()
    with filter_container:
        filter_columns = st.multiselect("Filter dataframe on", df.columns)

    for col in filter_columns:
        left, right = st.columns((1, 20))
        left.write("↳")

        if df[col].dtype == 'datetime64[ns]':
        #if len(df[col]) > 0 and isinstance(df[col][0], pd._libs.tslibs.timestamps.Timestamp):
            user_date_input = right.date_input(f"Values for {col}",
            value=(
                df[col].min(),
                df[col].max(),
                ),)
            
            if len(user_date_input) == 2:
                user_date_input = tuple(map(pd.to_datetime, user_date_input))
                start_date, end_date = user_date_input
                df = df.loc[df[col].between(start_date, end_date)]

        else:
            user_text_input = right.selectbox(
                f"{col}", df[col])
            #if user_text_input:
            df = df[df[col].str.contains(user_text_input)]

    for col in df.columns:
        if df[col].dtype == 'datetime64[ns]':
            df[col] = df[col].dt.strftime("%d-%m-%Y")

    return df

database_dir = os.getenv('CANCERMUTS_DATABASE')
if database_dir is None:
    database_dir = './database'

st.set_page_config(layout="wide",
    page_title="Cancermuts",
    page_icon="📖")

st.header("Welcome to Cancermuts!")
st.write("Navigate below to explore the data.")

try:
    show_table = pd.read_csv(os.path.join(database_dir, 'index_table.csv'))
except FileNotFoundError:
    st.write('No entries are currently available.')
    st.stop()

addFilter = st.checkbox("Add filter")

df = pd.DataFrame(show_table)

if addFilter:  
    df = filter_dataframe(df)

selection = get_selection(df)

if len(selection) != 1:
    invalid_selection = True
    protein=''
    data=''
else:
    invalid_selection = False
    protein = selection['Gene'].to_list()[0]
    data = os.path.join(database_dir, [filename for filename in os.listdir(database_dir) if filename.startswith(protein)][0])

st.download_button(label="Download dataset",
                disabled=invalid_selection,
                data=data,
                file_name=f'{protein}.csv',
                mime="text/csv",
                key='download-csv')

if invalid_selection == False:
        
    data = load_dataset(database_dir, protein)
    data = data.drop(['Unnamed: 0', 'Unnamed: 0.1'], axis=1)
    # start = 0
    start = data['aa_position'].iloc[0]
    # end = len(data)
    end = data['aa_position'].iloc[-1]

    def update_slider():
        st.session_state["slider"] = (st.session_state["nb_input1"], st.session_state["nb_input2"])

    slider_values = st.slider('Specify region of interest:', start, end, (start, end), key='slider')

    subcol1, subcol2 = st.columns(2)

    with subcol1:
        nb_input_val1 = st.number_input("Start", min_value=start, max_value=slider_values[1], value=slider_values[0], on_change=update_slider, key='nb_input1')

    with subcol2:
        nb_input_val2 = st.number_input("End", min_value=slider_values[0], max_value=end, value=slider_values[1], on_change=update_slider, key='nb_input2')

    start = nb_input_val1
    end = nb_input_val2

    start_i = min(data[data['aa_position'] == start].index.to_list())
    end_i = max(data[data['aa_position'] == end].index.to_list())
    region_specific_data = data[start_i:end_i+1].reset_index(drop=True)

    dataset, upset, revel, cancermuts = st.tabs(["Data", "UpSet plot", "REVEL distribution plot", "Cancermuts plot"])

    with dataset:
        st.write(f"Currently viewing: {protein}, in the AA range {start}-{end}")
        st.dataframe(region_specific_data)

    with upset:
        upset_data = region_specific_data[region_specific_data['sources'].notna()]
        upset_df = from_memberships(upset_data.sources.str.split(','), data=upset_data)
        upset_df.index.rename({'Manual annotations from mutations_clinvar.csv': 'clinvar',
                               'Manual annotations from clinvar.csv': 'clinvar'}, inplace=True)
        
        fig = plt.figure()
        plt.title(f'UpSet plot: {protein}, {start}-{end}')
        plt.axis('off')
        plot(upset_df, fig=fig)
        st.pyplot(fig=fig)

        download_plot_button(f'{protein}_upset_'+str(slider_values[0])+'-'+str(slider_values[1]))

    with revel:
        fig, ax = plt.subplots()
        hist_data = region_specific_data['REVEL_score'][region_specific_data['REVEL_score'].notna()].tolist()

        if len(hist_data) > 0:
            if type(hist_data[0]) == str:
                for i in range(len(hist_data)):
                    vals = hist_data[i].split(',')
                    if len(vals) > 1:
                        hist_data.pop(i)
                        for val in vals:
                            hist_data.append(val)
                hist_data = list(map(float, hist_data))

            nb_scores = len(hist_data)
            ax.hist(hist_data, bins=np.arange(0, 1, 0.05), ec='black')
            ax.set_xlabel("REVEL score")
            ax.set_ylabel("count")
            ax.set_title(f"REVEL score distribution plot: {protein} {start}-{end}")
            xticks = list(map(lambda x: round(x, 2), np.arange(0, 1, 0.05))) + [1.0]
            ax.set_xticks(xticks, labels=xticks, fontsize=7)
            ax.set_xlim(0,1)
            st.pyplot(fig)
            st.write(f"Total number of assigned REVEL scores in region {start}-{end}: {nb_scores}")

        else:
           st.write("No REVEL scores reported for this protein")

        download_plot_button(f'{protein}_revel_'+str(slider_values[0])+'-'+str(slider_values[1]))

    with cancermuts:
        vert_size = 10
        section = 50

        # for i in range(int(round(len(region_specific_data)/100))):
        #     vert_size += 2
        #     section -= 5

        # st.write(vert_size, section)

        if len(region_specific_data) <= 100:
            vert_size = 6
            section = 40
        elif len(region_specific_data) > 100 and len(region_specific_data) <= 200:
            vert_size = 10
            section = 40
        elif len(region_specific_data) > 200 and len(region_specific_data) <= 300:
            vert_size = 12
            section = 35
        elif len(region_specific_data) > 300 and len(region_specific_data) <= 400:
            vert_size = 17
            section = 30
        elif len(region_specific_data) > 400 and len(region_specific_data) <= 500:
            vert_size = 22
            section = 30
        elif len(region_specific_data) > 500 and len(region_specific_data) <= 600:
            vert_size = 30
            section = 25
        elif len(region_specific_data) > 600 and len(region_specific_data) <= 700:
            vert_size = 38
            section = 25

        tbl = Table()
        fig, ax = tbl.plot_metatable(region_specific_data, section_size=section, elm_y_ladder=(-0.2, -0.5, 5), rcParams={'font.size':8.0, 'font.sans-serif':['Arial']}, figsize=(10,vert_size))
        st.pyplot(fig=fig)

        download_plot_button(f'{protein}_cancermuts_'+str(slider_values[0])+'-'+str(slider_values[1]))
