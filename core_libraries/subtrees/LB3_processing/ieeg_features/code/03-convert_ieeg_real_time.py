# %%
# Imports
import numpy as np
import pandas as pd

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

tz = 'US/Eastern'
# path where lb3 data files are
root_path = '/mnt/local/gdrive/public/DATA/Human_Data/LB3_PIONEER/'
# %%
# Read metadata
metadata_path = "../data/ieeg_metadata.xlsx"
metadata = pd.read_excel(metadata_path)

# %%
for index, _ in metadata.iterrows():
    metadata.at[index, 'Start day/time'] = pd.Timestamp.combine(metadata['Start day'][index], metadata['Start time'][index]).tz_localize(tz='US/Eastern')
# %%
for index, row in metadata.iterrows():
    iEEG_times = row['Seizure events (iEEG time)']

    if np.isnan(iEEG_times):
        continue

    new_timestamp = pd.to_datetime(metadata['Start day/time'][index]) +  pd.Timedelta(iEEG_times, unit='s')
    metadata.at[index, 'Seizure events (real time)'] = new_timestamp

metadata.to_pickle("../data/ieeg_metadata_converted_sz_times.pkl")
metadata.to_csv("../data/ieeg_metadata_converted_sz_times.csv")
# %%

# %%
