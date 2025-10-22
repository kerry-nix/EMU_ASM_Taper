#%%
"""
This function provides a test case of pulling iEEG data
"""
# pylint: disable-msg=C0103
# pylint: disable-msg=C0301
#%%
# Imports
import pandas as pd
import numpy as np

import tools
import os

work_path = os.path.dirname(os.path.abspath(__file__))

# %%
with open(os.path.join(work_path, "config.json"), "rb") as f:
    config = pd.read_json(f, typ="series")

# iEEG_filename = "HUP172_phaseII"
# start_time_usec = 402580
# stop_time_usec = 402600
# electrodes = ["LE10", "LE11", "LH01", "LH02", "LH03", "LH04"]
iEEG_filename = "HUP212_phaseII"
start_time_usec = 402580
stop_time_usec = 402600
# %%
data, channels, fs = tools.get_ieeg_data(
    config.usr, config.pwd, iEEG_filename, start_time_usec, stop_time_usec,
)

# %%
# clean channel labels
clean_channels = tools.clean_labels(channels)
# data = tools.pre_whiten(data)
# find and return a boolean mask for non ieeg channels
non_ieeg_channels = tools.find_non_ieeg(clean_channels)

#%%
bad_chans, _ = tools.identify_bad_chan(data,fs)
#%%
keep_chans = (~non_ieeg_channels & ~bad_chans).squeeze()
data = data[:, keep_chans]
clean_channels = clean_channels[keep_chans]
#%% car
# data, car_labels = tools.common_average_ref(data, clean_channels)
# print(car_labels)
#%% bipolar
# data, bipolar_labels = tools.automatic_bipolar_montage(data, clean_channels)
# print(bipolar_labels)
#%% laplacian
locs = pd.read_csv(os.path.join(work_path,"tools","locs.csv"),header=None).values
locs = locs[keep_chans]
data, close_chs, laplacian_labels = tools.laplacian_reference(data,locs,5,clean_channels)
# print(close_chs)
print(data.shape)
print(laplacian_labels)
# %% Plot the data
# t_sec = np.linspace(start_time_usec, stop_time_usec, num=data.shape[0]) / 1e6
# fig, ax = tools.plot_iEEG_data(data, t_sec)
# fig.set_size_inches(18.5, 10.5)
# ax.set_title(iEEG_filename)
