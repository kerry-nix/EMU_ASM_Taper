# %%
# Imports
import os
from os.path import exists as ospe
from os.path import join as ospj

import pandas as pd
from ieeg.auth import Session

import warnings
warnings.filterwarnings("ignore")

username = 'pattnaik'
pwd_bin_path = "/mnt/leif/littlab/users/pattnaik/utils/pat_ieeglogin.bin"
with open(pwd_bin_path, "r") as f:
    s = Session(username, f.read())             


# path where lb3 data files are
root_path = "/gdrive/public/DATA/Human_Data/LB3_PIONEER"
# %%
# Read metadata
metadata_path = "/gdrive/public/DATA/Human_Data/LB3_PIONEER/patient_metadata.xlsx"
metadata = pd.read_excel(metadata_path)
# %%
hup_names_table = metadata[['LB3_id', 'iEEG_portal_names', 'Ignore']].drop_duplicates()

for _, row in hup_names_table.iterrows():
    if row['Ignore']:
        continue
    lb3_id = row['LB3_id']
    portal_name = row['iEEG_portal_names']
    save_path = ospj(root_path, lb3_id, "ieeg_features", f"{portal_name}.h5")
    if ospe(save_path):
        continue
    if not ospe(ospj(root_path, lb3_id, "ieeg_features")):
        os.makedirs(ospj(root_path, lb3_id, "ieeg_features"))
        
    print(portal_name)
    ds = s.open_dataset(portal_name)

    all_layers = ds.get_annotation_layers()

    layer_name = list(all_layers.keys())[0]
    expected_count = all_layers[layer_name]

    actual_count = 0
    max_results = None if expected_count < 100 else 100
    call_number = 0
    all_annotations_li = []
    while actual_count < expected_count:
        annotations = ds.get_annotations(
            layer_name, first_result=actual_count, max_results=max_results)
        call_number += 1
        actual_count += len(annotations)

        for annotation in annotations:
            first = pd.to_datetime(annotation.start_time_offset_usec, unit="us")
            duration = pd.to_datetime(annotation.end_time_offset_usec, unit="us") - first
            description = annotation.description
            all_annotations_li.append([first, duration, description])

        first = annotations[0].start_time_offset_usec
        last = annotations[-1].end_time_offset_usec
        description = annotations[0].description

        print("got", len(annotations), "annotations on call #",
                call_number, "covering", first, "usec to", last, "usec")

    all_annotations_df = pd.DataFrame(all_annotations_li)
    all_annotations_df.columns = ["Time", "Duration", "Description"]
    all_annotations_df.set_index("Time", inplace=True)

    all_annotations_df.to_hdf(save_path, key='annotations', mode='w')        

    print("Wrote all annotations to {}.h5".format(portal_name))
# %%
