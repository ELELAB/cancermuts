import os
import pandas as pd
from pathlib import Path
import streamlit_utils as su
from bioservices.uniprot import UniProt

datadir = "/data/raw_data/computational_data/cancermuts_data/"
entry_details_file = su.get_database_dir("example_entries.csv")
database_dir = su.get_database_dir("example_database")

entry_details_df = pd.read_csv(entry_details_file)

metatable = pd.DataFrame(columns=['Gene', 'Uniprot ID', 'Date', 'Source'])

for index, row in entry_details_df.iterrows():
    path = os.path.join(datadir, row['protein_class'], row['protein_name'], row['run_type'])
    #path = os.path.join(datadir, row['protein_class'], row['protein_name'].lower(), row['run_type'])

    if os.path.isdir(path):
        most_recent = pd.to_datetime("01-01-1900", format="%d-%m-%Y")
        for date in os.listdir(path):
            if os.path.isdir(os.path.join(path, date)):
                date_obj = pd.to_datetime(date, format="%d%m%Y")
                if date_obj > most_recent:
                    most_recent = date_obj
                
        dir_path = os.path.join(path, most_recent.strftime("%d%m%Y"))
        filename = [filename for filename in os.listdir(dir_path) if filename.startswith('metatable')][0]
        file_path = os.path.join(dir_path, filename)

        if os.path.isfile(file_path) and os.path.isfile(os.path.join(database_dir, filename)) != True:
            file = pd.read_csv(file_path)
            file.to_csv(os.path.join(database_dir, row['protein_name']))
            #file.to_csv(os.path.join(database_dir, row['protein_name'] + '_' + most_recent.strftime("%d%m%Y")))

            metatable_row = []
            metatable_row.append(row['protein_name'])
            u = UniProt(verbose=False)
            up_acc = u.search(row['protein_name']+"+and+taxonomy_id:9606", limit=3, columns="accession")
            metatable_row.append(up_acc.split('\n')[1])
            metatable_row.append(most_recent)
            metatable_row.append(row['Reference'])
            metatable.loc[len(metatable)] = metatable_row

    else:
        print(f"Entry not found: {path} is not a directory")

metatable.to_csv(os.path.join(database_dir, 'metatable.csv'), index=False)



    

