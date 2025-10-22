# %%
# Imports
from os.path import join as ospj

import pandas as pd
from ieeg.auth import Session

import warnings
warnings.filterwarnings("ignore")

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

username = 'pattnaik'
pwd_bin_path = "/mnt/leif/littlab/users/pattnaik/utils/pat_ieeglogin.bin"
with open(pwd_bin_path, "r") as f:
    s = Session(username, f.read())             

keywords = ["seizure", "sz", "ictal", "event", "jerks", "auras", "episodes", "myoclon", "spell", "EEC", "UEO", "offset"]
keywords_regex = "|".join(keywords)
# for making it case insenstive
keywords_regex = "^(?i)" + keywords_regex

# path where lb3 data files are
root_path = "/gdrive/public/DATA/Human_Data/LB3_PIONEER"
# %%
# Read metadata
metadata_path = "/gdrive/public/DATA/Human_Data/LB3_PIONEER/patient_metadata.xlsx"
metadata = pd.read_excel(metadata_path)
# %%
hup_names_table = metadata[['LB3_id', 'iEEG_portal_names']].drop_duplicates()

for _, row in hup_names_table.iterrows():
    lb3_id = row['LB3_id']
    portal_name = row['iEEG_portal_names']

    fname = ospj(root_path, lb3_id, "ieeg_features", f"{portal_name}.h5")
    key = 'annotations'

    all_annotations = pd.read_hdf(fname, key)

    mask = all_annotations['Description'].str.contains(keywords_regex)
    mask[mask.isna()] = False

    filtered_annotations = all_annotations[mask]

    filter_pct = 1 - len(filtered_annotations) / len(all_annotations)
    print("Filtered {:10.2f}% of all annotations".format(filter_pct*10))

    t_sec = [i.timestamp() for i in filtered_annotations.index]
    filtered_annotations['IEEG Time'] = t_sec
    try:
        filtered_annotations.to_excel(ospj(root_path, f"{lb3_id}/ieeg_features/{portal_name}_filtered_annotations.xlsx"))
        filtered_annotations.to_hdf(ospj(root_path, lb3_id, "ieeg_features", f"{portal_name}.h5"), key='annotations_filtered', mode='a')        

        print("Wrote filtered annotations to {}.h5".format(portal_name))
    except:
        print("Did not write filtered annotations to {}.h5".format(portal_name))
        continue

# %%
for _, row in hup_names_table.iterrows():
    lb3_id = row['LB3_id']
    portal_name = row['iEEG_portal_names']

    if lb3_id != 'LB3_002_phaseI':
        continue

    fname = ospj(root_path, lb3_id, "ieeg_features", f"{portal_name}.h5")

    filtered_annotations = pd.read_hdf(fname, 'annotations_filtered')
    filtered_annotations.to_excel(f"../../../{lb3_id}/ieeg_features/{portal_name}_filtered_annotations.xlsx")
    
# %%
