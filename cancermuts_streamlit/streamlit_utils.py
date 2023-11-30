# streamlit_utils.py for cancermuts
# (c) 2023 Alberte Heering Estad <ahestad@outlook.com>
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

import streamlit as st
import os
import pandas as pd

def get_selection(df: pd.DataFrame) -> pd.DataFrame:
    df_with_selections = df.copy()
    df_with_selections.insert(0, "Select", False)

    edited_df = st.data_editor(
        df_with_selections,
        hide_index=True,
        column_config={"Select": st.column_config.CheckboxColumn(required=True)},
        disabled=df.columns,
    )
    # Filter the dataframe using the temporary column, then drop the column
    selected_rows = edited_df[edited_df.Select]
    return selected_rows.drop('Select', axis=1)

def get_database_dir(var_name='cancermuts_data', default_dir='./database'):
    dir_name = os.getenv(var_name)
    if dir_name is None:
        return default_dir
    return var_name

# def get_database_dir(database):
#     return os.path.join(os.getcwd(), database)

def load_dataset(data_dir, protein):
    #return pd.read_csv(os.path.join(data_dir, f'{protein}_2*'))
    return pd.read_csv(os.path.join(data_dir, [filename for filename in os.listdir(data_dir) if filename.startswith(protein)][0]))

def load_main_table(data_dir):
    return pd.read_csv(os.path.join(data_dir, 'index_table.csv'))