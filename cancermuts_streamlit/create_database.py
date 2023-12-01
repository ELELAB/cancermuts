# create_database.py for cancermuts
# (c) 2023 Alberte Heering Estad <ahestad@outlook.com>, Matteo Tiberti
# Danish Cancer Institute
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

import os
import pandas as pd
import streamlit_utils as su
from bioservices.uniprot import UniProt
import collections
import argparse

argparser = argparse.ArgumentParser(description='Create a database of cancermuts entries for the web application')

argparser.add_argument('-d', '--datadir', type=str, help='Path to the directory containing cancermuts data', required=True)
argparser.add_argument('-e', '--entry_details_file', type=str, help='Path to the file containing the entry details', required=True)
argparser.add_argument('-o', '--output_dir', type=str, help='Path to the output directory', required=True)

args = argparser.parse_args()

u = UniProt(verbose=False)

datadir = args.datadir
entry_details_file = args.entry_details_file
database_dir = args.output_dir

entry_details_df = pd.read_csv(entry_details_file)

rows = collections.defaultdict(list)

for index, row in entry_details_df.iterrows():
    path = os.path.join(datadir, row['protein_class'], row['protein_name'], row['run_type'])
    #path = os.path.join(datadir, row['protein_class'], row['protein_name'].lower(), row['run_type'])

    if not os.path.isdir(path):
        print(f"Entry not found: {path} is not a directory")
        continue

    most_recent = pd.to_datetime("01-01-1900", format="%d-%m-%Y")
    for date in os.listdir(path):
        if os.path.isdir(os.path.join(path, date)):
            date_obj = pd.to_datetime(date, format="%d%m%Y")
            if date_obj > most_recent:
                most_recent = date_obj

    dir_path = os.path.join(path, most_recent.strftime("%d%m%Y"))
    filename = [filename for filename in os.listdir(dir_path) if filename.startswith('metatable') and filename.endswith('.csv')]

    if len(filename) > 1:
        print("More than one Cancermuts metatable was found in the path: " + dir_path)
        continue

    file_path = os.path.join(dir_path, filename[0])

    if os.path.isfile(file_path) and not os.path.isfile(os.path.join(database_dir, filename[0])):
        file = pd.read_csv(file_path)
        file.to_csv(os.path.join(database_dir, row['protein_name']))

        up_acc = u.search(row['protein_name']+"+and+taxonomy_id:9606", limit=3, columns="accession")

        rows['Gene'].append(row['protein_name'])
        rows['Uniprot ID'].append(up_acc.split('\n')[1])
        rows['Date'].append(most_recent)
        rows['Source'].append(row['Reference'])

df = pd.DataFrame.from_dict(rows)
df.to_csv(os.path.join(database_dir, 'index_table.csv'), index=False)





