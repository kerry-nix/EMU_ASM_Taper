# %%
# Imports
# %load_ext autoreload
# %autoreload 2
import numpy as np
import pandas as pd
from tools import get_bandpower
import os
from ieeg.auth import Session
import pickle
import matplotlib.pyplot as plt
from datetime import datetime, date
import scipy.cluster.hierarchy as shc
from sklearn.cluster import KMeans
from scipy.stats import kruskal
import scikit_posthocs as sp

def filter_df_by_times(df, t1, t2):
    ind1 = df.index.searchsorted(t1)
    ind2 = df.index.searchsorted(t2)
    print(ind1, ind2)
    return(df.iloc[ind1:ind2-1, :])

pd.reset_option('display.max_columns')
pd.reset_option('display.max_rows')

PLOT_FLAG = True

tz = 'US/Eastern'
time_before_after_sec = int(2.5 * 60)
posix = pd.Timestamp("1970-01-01", tz=tz)
# path where lb3 data files are
root_path = '/mnt/local/gdrive/public/DATA/Human_Data/LB3_PIONEER/'
# %%
# Read metadata
metadata_path = "../data/ieeg_metadata_converted_sz_times.xlsx"
metadata = pd.read_excel(metadata_path)
metadata_path = "../data/ieeg_metadata_converted_sz_times.pkl"
metadata = pd.read_pickle(metadata_path)

username = 'pattnaik'
pwd_bin_path = "/gdrive/public/USERS/pattnaik/utils/pat_ieeglogin.bin"
with open(pwd_bin_path, "r") as f:
    s = Session(username, f.read())             # start an IEEG session with your username and password. TODO: where should people put their ieeg_pwd.bin file?

# # %%
# for _, row in metadata.iterrows():
#     if np.isnan(row['Seizure events (iEEG time)']):
#         continue

#     if row['LB3_id'] == 'LB3_001_phaseII':
#         continue

#     lb3_id = row['LB3_id']
#     portal_name = row['iEEG_portal_names']
        
#     print(lb3_id, row['Seizure events (iEEG time)'])

#     posix_circadian_timeshift = pd.Timedelta(datetime.combine(date.min, row['Start time']) - datetime.min)
    
#     sz_data_file = "../../../{}/{}_seizure_data.pkl".format(lb3_id, lb3_id)
#     with open(sz_data_file, 'rb') as f:
#         seizure_data = pickle.load(f)
#     # assert(seizure_data['LB3_id'] == lb3_id)

#     for sz in seizure_data['seizure_data']:
#         sz_timestamps = [(i['seizure_timestamp'] - row['Start day/time']).total_seconds() for i in seizure_data['seizure_data']]
#         sz_index = sz_timestamps.index(row['Seizure events (iEEG time)'])
#         # acc_data = sz['wearables_around_surveys'].rename(columns={"accelerometerAccelerationX(G)": "acc_x", "accelerometerAccelerationY(G)": "acc_y", "accelerometerAccelerationZ(G)": "acc_z"})
#         # acc_data = acc_data.set_index('timestamp')
#         # # get and format accelerometer data to posix time with time of day preserved
#         # acc_data = acc_data[acc_data['acc_x'].notna()]
#         # posix_timeshift = row['Start day/time'] - posix
#         # indices = ([i - posix_timeshift + posix_circadian_timeshift for i in acc_data.index])
#         # acc_data.index = indices

#     sz = seizure_data['seizure_data'][sz_index]
#     # get timestamps
#     seizure_timestamp = pd.Timestamp(sz['seizure_timestamp'])
#     before_survey_timestamp = sz['survey_before_seizure']['primary_survey_timestamp']
#     after_survey_timestamp = sz['survey_after_seizure']['primary_survey_timestamp']

#     # convert timestamp to iEEG time
#     seizure_ieeg = (seizure_timestamp - row['Start day/time']).total_seconds()
#     before_survey_ieeg = (before_survey_timestamp - row['Start day/time']).total_seconds()
#     after_survey_ieeg = (after_survey_timestamp - row['Start day/time']).total_seconds()

#     # convert iEEG time to posix time
#     seizure_posix = posix + pd.Timedelta(seizure_ieeg, unit='s')
#     before_survey_posix = posix + pd.Timedelta(before_survey_ieeg, unit='s')
#     after_survey_posix = posix + pd.Timedelta(after_survey_ieeg, unit='s')

#     # convert posix time to posix time with time of day preserved
#     seizure_posix_circadian = seizure_posix + posix_circadian_timeshift
#     before_survey_circadian = before_survey_posix + posix_circadian_timeshift
#     after_survey_circadian = after_survey_posix + posix_circadian_timeshift

    
#     sz_bandpowers = get_bandpower(
#         username, 
#         pwd_bin_path, 
#         portal_name, 
#         seizure_ieeg - time_before_after_sec, 
#         seizure_ieeg + row['Seizure duration'] + time_before_after_sec,
#         bands_kw='broad'
#         )

#     before_bandpowers = get_bandpower(
#         username, 
#         pwd_bin_path, 
#         portal_name, 
#         before_survey_ieeg - time_before_after_sec, 
#         before_survey_ieeg  + time_before_after_sec,
#         bands_kw='broad'
#         )
#     after_bandpowers = get_bandpower(
#         username, 
#         pwd_bin_path, 
#         portal_name, 
#         after_survey_ieeg - time_before_after_sec, 
#         after_survey_ieeg  + time_before_after_sec,
#         bands_kw='broad'
#         )

#     # shift index to posix time with time of day preserved
#     before_bandpowers.index = before_bandpowers.index.tz_localize(tz=tz).shift(freq=posix_circadian_timeshift)
#     sz_bandpowers.index = sz_bandpowers.index.tz_localize(tz=tz).shift(freq=posix_circadian_timeshift)
#     after_bandpowers.index = after_bandpowers.index.tz_localize(tz=tz).shift(freq=posix_circadian_timeshift)

#     foldername = "sz_{:.0f}".format(row['Seizure events (iEEG time)'])
#     fpath = '../../../{}/ieeg_features'.format(row['LB3_id'])
#     if not os.path.exists(os.path.join(fpath, foldername)):
#         os.makedirs(os.path.join(fpath, foldername))

#     before_bandpowers.to_pickle(os.path.join(fpath, foldername, "before_clip_bandpowers.pkl"))
#     sz_bandpowers.to_pickle(os.path.join(fpath, foldername, "seizure_clip_bandpowers.pkl"))
#     after_bandpowers.to_pickle(os.path.join(fpath, foldername, "after_clip_bandpowers.pkl"))

# %%
# If returning, start here
for _, row in metadata.iterrows():
    if np.isnan(row['Seizure events (iEEG time)']):
        continue
    if row['LB3_id'] == 'LB3_001_phaseII':
        continue

    # Extract variables from metadata
    lb3_id = row['LB3_id']
    start_time = row['Start time']
    iEEG_time = row['Seizure events (iEEG time)']

    print(lb3_id, row['Seizure events (iEEG time)'])

    # Set circadian timeshift for alignment
    posix_circadian_timeshift = pd.Timedelta(datetime.combine(date.min, start_time) - datetime.min)

    # Load in datafile which contains wearable and survey (generated by K. Xie and J. Panchal)
    sz_data_file = "../../../{}/{}_seizure_data.pkl".format(lb3_id, lb3_id)
    with open(sz_data_file, 'rb') as f:
        seizure_data = pickle.load(f)
    # assert(seizure_data['LB3_id'] == lb3_id)

    # Set paths and make sure directories are made
    foldername = "sz_{:.0f}".format(iEEG_time)
    fpath = '../../../{}/ieeg_features'.format(lb3_id)

    figure_path = '../../../{}/figures'.format(lb3_id)
    if not os.path.exists(os.path.join(figure_path, foldername)):
        os.makedirs(os.path.join(figure_path, foldername))

    # Prepare data in lists to iterate and create figures
    before_bandpowers = pd.read_pickle(os.path.join(fpath, foldername, "before_clip_bandpowers.pkl"))
    sz_bandpowers = pd.read_pickle(os.path.join(fpath, foldername, "seizure_clip_bandpowers.pkl"))
    after_bandpowers = pd.read_pickle(os.path.join(fpath, foldername, "after_clip_bandpowers.pkl"))

    all_bandpowers = [before_bandpowers, sz_bandpowers, after_bandpowers]
    wearables_keys = ['wearables_around_survey_before_seizure', 'wearables_around_seizure', 'wearables_around_survey_after_seizure']

    periods = ['before', 'seizure', 'after']

    vlines = [
        seizure_data['seizure_data'][0]['survey_before_seizure']['primary_survey_timestamp'],
        None,
        seizure_data['seizure_data'][0]['survey_after_seizure']['primary_survey_timestamp']
        ]

    # Loop to create figures for before, during, and after seizure
    for period, wearables_key, bandpowers, vline in zip(periods, wearables_keys, all_bandpowers, vlines):
        # Reset wearables index to aligned timestamps
        patient_seizure_times = [(i['seizure_timestamp'] - row['Start day/time']).total_seconds() for i in seizure_data['seizure_data']]
        seizure_idx = patient_seizure_times.index(iEEG_time)
        wearables_data = seizure_data['seizure_data'][seizure_idx][wearables_key]
        # wearables_data = wearables_data.rename(columns={"accX": "acc_x", "accelerometerAccelerationY(G)": "acc_y", "accelerometerAccelerationZ(G)": "acc_z"})
        wearables_data = wearables_data.set_index('timeStamp')
        wearables_data.index = [i - row['Start day/time'] + posix + posix_circadian_timeshift for i in wearables_data.index]
        
        n_clusters = 4
        kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit_predict(bandpowers)
        kmeans_df = pd.DataFrame(kmeans, index=bandpowers.index, columns=['State'])
        
        # Find the average energy in each state
        states_avg_energy = []
        past_index = posix
        for ind, kmeansrow in kmeans_df.iterrows():
            current_index = ind

            wearables_data_clip = wearables_data[(wearables_data.index > past_index) & (wearables_data.index < current_index)]

            # calculate energy
            wearables_data_clip['Average Energy'] = np.sqrt(wearables_data_clip['accX']*wearables_data_clip['accX'] + wearables_data_clip['accY']*wearables_data_clip['accY'] + wearables_data_clip['accZ']*wearables_data_clip['accZ'])

            avg_energy = wearables_data_clip['Average Energy'].mean()
            states_avg_energy.append((kmeansrow['State'], avg_energy))

            past_index = current_index
        states_avg_energy = pd.DataFrame(states_avg_energy, columns = ["State", "Average Energy"])        
        kw_test = kruskal(*[group["Average Energy"].values for name, group in states_avg_energy.groupby("State")])

        ph_dunn = sp.posthoc_dunn([group["Average Energy"].values for name, group in states_avg_energy.groupby("State")], p_adjust = 'bonferroni')

        if PLOT_FLAG:
            # if vline:
            #     vline = vline - row['Start day/time'] + posix + posix_circadian_timeshift

            # Visualize dendrogram for clusters
            ####
            fig, ax = plt.subplots(figsize=(10, 7))
            ax.set_title(f"Dendrograms for {period}, {lb3_id} {row['Seizure events (iEEG time)']}")  
            dend = shc.dendrogram(shc.linkage(bandpowers, method='ward'))
            ax.set_ylabel("Linkage distance")
            plt.savefig(os.path.join(figure_path, foldername, "dendrogram_period-{}.png".format(period)))
            plt.savefig(os.path.join(figure_path, foldername, "dendrogram_period-{}.svg".format(period)))
            plt.close()
            ####

            # Plot average bandpower across all channels aligned with accelerometer
            ####
            fig, axes = plt.subplots(nrows = 2, figsize=(16, 8), dpi=150, sharex=True)
            axes[0].plot(wearables_data)
            axes[0].legend(['x', 'y', 'z'])
            axes[0].set_ylabel("acc (mm^3)")
            axes[1].plot(np.mean(bandpowers, axis=1), color='k')
            axes[1].set_title("Average bandpower across all channels")
            axes[1].set_ylabel("Relative bandpower (out of 1)")        
            axes[1].set_xlabel("Time (time of day preserved)")

            fig.suptitle(f"{lb3_id} {row['Seizure events (iEEG time)']} {period}")

            # if vline:
            #     axes[0].axvline(vline, color='k', ls='--')
            #     axes[1].axvline(vline, color='k', ls='--')
            plt.savefig(os.path.join(figure_path, foldername, "average_bandpower_period-{}.png".format(period)))
            plt.savefig(os.path.join(figure_path, foldername, "average_bandpower_period-{}.svg".format(period)))
            plt.close()
            ####

            # Cluster each time period and plot states aligned with accelerometer
            ####
            fig, axes = plt.subplots(nrows = 2, figsize=(16, 8), dpi=150, sharex=True)
            axes[0].plot(wearables_data)
            axes[0].legend(['x', 'y', 'z'])
            axes[0].set_ylabel("acc (m/s^2)")
            axes[1].plot(pd.DataFrame(kmeans, index=bandpowers.index))
            axes[1].set_yticks(np.arange(n_clusters))
            axes[1].set_yticklabels(np.arange(n_clusters))
            axes[1].set_title("Clustering State (n = {})".format(n_clusters))
            axes[1].set_ylabel("State (a.u.)")
            axes[1].set_xlabel("Time (time of day preserved)")
            fig.suptitle(f"{lb3_id} {row['Seizure events (iEEG time)']} {period}")
            # if vline:
            #     axes[0].axvline(vline, color='k', ls='--')
            #     axes[1].axvline(vline, color='k', ls='--')

            plt.savefig(os.path.join(figure_path, foldername, "kmeans_clusters-{}_period-{}.png".format(n_clusters, period)))
            plt.savefig(os.path.join(figure_path, foldername, "kmeans_clusters-{}_period-{}.svg".format(n_clusters, period)))
            plt.close()


            fig, ax = plt.subplots()
            states_avg_energy.boxplot(column='Average Energy', by='State', ax=ax)
            ax.text(0.9, 1.1, f"Kruskal-Wallis = {kw_test[0]:.2f}\np = {kw_test[1]:.2E}", horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes)
            fig.suptitle("")
            ax.set_title(f"{lb3_id} {row['Seizure events (iEEG time)']} {period}")
            ax.set_ylabel("Average Energy")

            plt.savefig(os.path.join(figure_path, foldername, "states_avg_energy-{}_period-{}.png".format(n_clusters, period)), bbox_inches='tight')
            plt.close()
# %%

# %%
