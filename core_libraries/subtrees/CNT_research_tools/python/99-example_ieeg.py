#%%
'''
This function provides a test case of pulling iEEG data
'''
# pylint: disable-msg=C0103
# pylint: disable-msg=C0301
#%%
# Imports
import pandas as pd
import numpy as np

import tools

# %%
with open("../config.json", 'rb') as f:
    config = pd.read_json(f, typ='series')

iEEG_filename = "HUP172_phaseII"
start_time_usec = 402580 * 1e6
stop_time_usec = 402600 * 1e6
electrodes = ["LE10","LE11","LH01","LH02","LH03","LH04"]

# %%
data, fs = tools.get_iEEG_data(config.usr, config.pwd, iEEG_filename, start_time_usec, stop_time_usec)

# %%
print(data.columns)

# %%
# clean channel labels
clean_channels = tools.clean_labels(data.columns)
# find and return a boolean mask for non ieeg channels
non_ieeg_channels = tools.find_non_ieeg(clean_channels)

data.columns = clean_channels
data = data.iloc[:, ~non_ieeg_channels]

#%%
print(data.columns)

# %% Plot the data
t_sec = np.linspace(start_time_usec, stop_time_usec, num=data.shape[0]) / 1e6
fig, ax = tools.plot_iEEG_data(data, t_sec)
fig.set_size_inches(18.5, 10.5)
ax.set_title(iEEG_filename)
fig.show()
