# %%
# Imports
import numpy as np
import multiprocessing # todo
from os.path import join as ospj
import sys
import pandas as pd
from ieeg.auth import Session

from numpy.lib.stride_tricks import sliding_window_view
from scipy.signal import iirnotch, filtfilt, butter

import warnings
warnings.filterwarnings("ignore")

# username = 'pattnaik'
# pwd_bin_path = "/gdrive/public/USERS/pattnaik/utils/pat_ieeglogin.bin"
# with open(pwd_bin_path, "r") as f:
#     s = Session(username, f.read())             # start an IEEG session with your username and password. TODO: where should people put their ieeg_pwd.bin file?


# # path where lb3 data files are
# root_path = '/mnt/local/gdrive/public/DATA/Human_Data/LB3_PIONEER/'
# # %%
# # Read metadata
# metadata_path = "../data/ieeg_metadata.xlsx"
# metadata = pd.read_excel(metadata_path)

# %%
# Bands
delta_band = [1, 4]
theta_band = [4, 8]
alpha_band = [8, 13]
beta_band = [13, 30]
gamma_band = [30, 70]
high_gamma_band = [70, 120]
broad_band = [5, 115]

bands = [
    delta_band,
    theta_band,
    alpha_band,
    beta_band,
    gamma_band,
    high_gamma_band,
    broad_band
]

band_names = [
    "delta",
    "theta",
    "alpha",
    "beta",
    "gamma",
    "high_gamma",
    "broad"
]

# %%
# Parameters
# size of window for each downloaded data chunk
# data_pull_min = 5
# # window size for each bandpower calculation (seconds)
# win_size_sec = 5
# # window stride (how many seconds overlap) for each bandpower calculation (seconds)
# win_stride_sec = 2.5

# %%
# for _, row in metadata.iterrows():
#     lb3_id = row['LB3_id']
#     iEEG_portal_names = [x.strip() for x in row['iEEG_portal_names'].split(',')]
    
    
#     for portal_name in iEEG_portal_names:


def get_bandpower(username, pwd_bin_path, portal_name, start_time_sec, end_time_sec, bands_kw='all', channels='all', data_pull_min=5, win_size_sec=5, win_stride_sec=2.5):
    from tqdm import tqdm
    from tools import get_iEEG_data, bandpower
    with open(pwd_bin_path, "r") as f:
        s = Session(username, f.read())
    ds = s.open_dataset(portal_name)

    clip_duration_sec =  end_time_sec - start_time_sec
    clip_duration_min = clip_duration_sec / 60

    # how many 5 minute data pulls are there?
    n_iter = int(np.ceil(clip_duration_min / data_pull_min))

    if channels == 'all':
        channel_labels = ds.get_channel_labels()

    if bands_kw == 'broad':
        bands = [[5, 115]]
        band_names = ['broad']
    elif bands_kw == 'all':
        delta_band = [1, 4]
        theta_band = [4, 8]
        alpha_band = [8, 13]
        beta_band = [13, 30]
        gamma_band = [30, 70]
        high_gamma_band = [70, 120]
        broad_band = [5, 115]

        bands = [
            delta_band,
            theta_band,
            alpha_band,
            beta_band,
            gamma_band,
            high_gamma_band,
            broad_band
        ]

        band_names = [
            "delta",
            "theta",
            "alpha",
            "beta",
            "gamma",
            "high_gamma",
            "broad"
        ]

    elif bands_kw != 'all':
        sys.exit('invalid band argument, should be broad or all')
    band_powers = None
    for i in tqdm(range(n_iter)):
        # This is for pulling ALL data, would take up too much memory
        # start_usec = i * (data_pull_min * 60 * 1e6)
        # end_usec = (i + 1) * (data_pull_min * 60 * 1e6)

        start_usec = start_time_sec * 1e6 + i * (5 * 60 * 1e6)
        if i == n_iter - 1:
            end_usec = end_time_sec * 1e6
        else:
            end_usec = start_usec + (i + 1) * (5 * 60 * 1e6)

        data, fs = get_iEEG_data(username, pwd_bin_path, portal_name, start_usec, end_usec, select_electrodes=channel_labels)

        time = np.linspace(start_usec, end_usec, len(data), endpoint=False)
        data.index = pd.to_datetime(time, unit='us')

        # extract dims
        n_samples = np.size(data, axis=0)
        n_channels = np.size(data, axis=1)

        assert(n_channels == len(channel_labels))

        win_size_ind = int(win_size_sec * fs)
        win_stride_ind = int(win_stride_sec * fs)

        win_size_ind = int(win_size_sec * fs)
        win_stride_ind = int(win_stride_sec * fs)

        sl_win = sliding_window_view(np.arange(n_samples), win_size_ind)[::win_stride_ind, :]

        # nan check
        nan_mask = np.ones(n_samples, dtype=bool)
        for win_inds in sl_win:
            if np.sum(np.isnan(data.iloc[win_inds, :]), axis=0).any():
                nan_mask[win_inds] = False
        signal_nan = data[nan_mask]

        # extract dims, again since nans may have reduced window size
        if np.sum(~nan_mask) > 0:
            n_samples = np.size(signal_nan, axis=0)
            n_channels = np.size(signal_nan, axis=1)
            sl_win = sliding_window_view(np.arange(n_samples), win_size_ind)[::win_stride_ind, :]

        # artifact rejection
        nan_mask = np.ones(n_samples, dtype=bool)
        for win_inds in sl_win:
            if (np.sum(np.abs(signal_nan.iloc[win_inds, :]), axis=0) < 1/12).any():
                nan_mask[win_inds] = False
            if (np.sqrt(np.sum(np.diff(signal_nan.iloc[win_inds, :]))) > 15000).any():
                nan_mask[win_inds] = False
            
        signal_minus_artifact = signal_nan[nan_mask]

        if len(signal_minus_artifact) == 0:
            continue
        
        # remove 60Hz noise
        f0 = 60.0  # Frequency to be removed from signal (Hz)
        Q = 30.0  # Quality factor
        b, a = iirnotch(f0, Q, fs)
        signal_filt = filtfilt(b, a, signal_minus_artifact, axis=0)

        # bandpass between 1 and 120Hz
        bandpass_b, bandpass_a = butter(3, [1, 120], btype='bandpass', fs=fs)
        signal_filt = filtfilt(bandpass_b, bandpass_a, signal_filt, axis=0)

        # format resulting data into pandas DataFrame
        signal_filt = pd.DataFrame(signal_filt, columns=signal_nan.columns, index=signal_nan.index)
        
        # re-reference the signals using common average referencing
        signal_ref = signal_filt.subtract(signal_filt.mean(axis=1), axis=0)

        all_bandpowers = None
        for name, band in zip(band_names, bands):
            sl_win = sliding_window_view(np.arange(signal_ref.shape[0]), win_size_ind)[::win_stride_ind, :]
            n_windows = sl_win.shape[0]

            # calculate bandpower
            power_mat = np.zeros((n_windows, n_channels))
            for ind, win_inds in enumerate(sl_win):
                power_mat[ind, :] = bandpower(signal_ref.iloc[win_inds, :], fs, band, 1, relative=True)

            if all_bandpowers is None:
                all_bandpowers = power_mat
            else:
                all_bandpowers = np.hstack((all_bandpowers, power_mat))
        

        col_names = [[ "{} {}".format(elec, band) for band in band_names for elec in signal_ref.columns ]]
        indices = signal_ref.index[sl_win[:, -1]]
        bandpower_df = pd.DataFrame(all_bandpowers, index=indices, columns=col_names)

        if band_powers is None:
            band_powers = bandpower_df
        else:
            band_powers = pd.concat([band_powers, bandpower_df], axis=0)

    return band_powers
    # band_powers.to_hdf(ospj(root_path, lb3_id, "ieeg_features", "{}.h5".format(portal_name)), key='bandpower', mode='w')
# %%
