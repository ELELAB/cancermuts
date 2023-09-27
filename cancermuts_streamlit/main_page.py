import streamlit as st
import os
import pandas as pd
from streamlit_utils import *
import matplotlib.pyplot as plt
from upsetplot import plot
from upsetplot import from_memberships
import numpy as np

def add_padding(amount: int):
    for i in range(amount):
        st.write('')

database_dir = get_database_dir('example_database')

st.set_page_config(layout="wide",
    page_title="Cancermuts",
    page_icon="📖")

st.header("Welcome to Cancermuts!")
st.write("Navigate below to view and download datasets.")

try:
    show_table = load_main_table(database_dir)
except FileNotFoundError:
    st.write('No entries are currently available.')
    st.stop()

df = pd.DataFrame(show_table)

selection = get_selection(filter_dataframe(df))

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
    start = 0
    end = len(data)

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

    region_specific_data = data[start:end+1]

    dataset, upset, revel = st.tabs(["Data", "UpSet plot", "REVEL distribution plot"])

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


