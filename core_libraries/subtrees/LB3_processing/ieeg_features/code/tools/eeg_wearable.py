import tools
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
from scipy.signal import iirnotch, filtfilt, butter
import importlib
import sys
import os
sys.path.insert(1, os.path.join("./tools/"))

# reload library
importlib.reload(sys.modules['helpers'])
from helpers import *

from tools import get_iEEG_data
from get_iEEG_data import *

def eeg_wearable(iEEG_filename, start_time, duration, subject_name, metadata_subjectname, channel_num1, channel_num2, channel_num3, channel_num4):

    # metadata_subjectname is often the same as subject_name
        # only not true when recording was interrupted, resulting in: subject_name_D01, subject_name_D02, etc... 
        # for metadata_subjectname 

    # start time is the listed staryt time of the seizre: the function will later subtract 60 from this in order to start showing it from a minute before the onset of the seizure

    # channel numbers to be included in the final plot MUST be indicated: look at iEEG.org or the first plot from obtaining  
    # the iEEG data in order to determine which channels show the most seiaure acyivity for detection 
            # --> *** Write the channel number as it is in iEEG... the function already changes the number to be indexed in Python syntax (starting from 0)

    #########################################
    
    ### ### FIRST, obtaining iEEG sample

    # tools & libraries to be used


    # # #

    # for removing 60 Hz noise later on...
    F0 = 60.0  # Frequency to be removed from signal (Hz)
    Q = 30.0  # Quality factor


    # bipolar reference function for obtaining bipolar data (iEEG sample)

    def _bipolar_reference(data_arg):
        (n_time_samples, n_channels) = data_arg.shape
        # separate contact names
        leads = []
        contacts = []
        for i in data_arg.columns:
            if i in ['C3', 'C4', 'ECG1', 'ECG2']:
                n_channels = n_channels - 1
                continue

            M = re.match(r"(\D+)(\d+)", i)
            if M is None:
                n_channels = n_channels - 1
                continue
            leads.append(M.group(1).replace("EEG", "").strip())
            contacts.append(int(M.group(2)))

        leads_contacts = [f"{i}{j:02d}" for i, j in zip(leads, contacts)]
        col_names = []
        data_entries = []
        # find montage channel before
        for index in range(n_channels - 1):
            lead = leads[index]
            contact = contacts[index]

            next_lead_contact = f"{lead}{(contact + 1):02d}"

            try:
                next_index = leads_contacts.index(next_lead_contact)
            except ValueError:
                continue

            col_names.append(f"{leads_contacts[index]}-{leads_contacts[next_index]}")
            data_entries.append(data_arg.iloc[:, index] - data_arg.iloc[:, next_index])
        
        data_entries = np.array(data_entries).T
        bipolar_data = pd.DataFrame(data_entries, columns=col_names, index=data_arg.index)

        return bipolar_data


    # iEEG login
    username = 'igvilla'
    pwd_bin_file = '/gdrive/public/USERS/igvilla/igv_ieeglogin.bin'


    # cutting time frame:
    start_time_sec = start_time - 60    # showing starting 1 minute before the seizure
    end_time_sec = start_time + duration

    # accessing iEEG file name & specific time clip:
    iEEG_filename = iEEG_filename
    start_time_usec = start_time_sec * 1e6
    end_time_usec = end_time_sec * 1e6

    # obtaining iEEG data
    data, fs = get_iEEG_data(username, pwd_bin_file, iEEG_filename, start_time_usec, end_time_usec, select_electrodes=np.arange(25))

    # extract dimensions
    (n_samples, n_channels) = data.shape

    # set time array
    t_sec = np.linspace(start_time_sec, end_time_sec, n_samples)

    # remove 60Hz noise
    b, a = iirnotch(F0, Q, fs)
    signal_filt = filtfilt(b, a, data, axis=0)

    # bandpass between 1 and 120Hz
    bandpass_b, bandpass_a = butter(3, [1, 120], btype='bandpass', fs=fs)
    signal_filt = filtfilt(bandpass_b, bandpass_a, signal_filt, axis=0)

    # format resulting data into pandas DataFrame
    signal_filt = pd.DataFrame(signal_filt, columns=data.columns)
    signal_filt.index = pd.to_timedelta(t_sec, unit="S")

    # defing the variable to be global, so that I can access this variable in the plot function
 #   global signal_ref

    signal_ref = _bipolar_reference(signal_filt)

    ### return signal_ref ###

    # create an hdf file out of the EEG data 
    #signal_ref.to_hdf(f"../data/eeg_pt-{iEEG_filename}_start-{start_time_sec}_end-{end_time_sec}.h5", key='ieeg')

    # plotting the EEG channels (ONLY) to observe which channels show more seizure activity 
    plt.figure(figsize=(10, 100)); plt.plot(signal_ref.iloc[:, 0:signal_ref.shape[1]] + np.arange(signal_ref.shape[1])*-1000)

    #########################################

    ### ### SECOND: overall, to obtain the wearables data & EEG data from the step above to plot


    # obtaining subject name & data
    subject_name = subject_name
    metadata = pd.DataFrame(json.load(open("./tools/subject_metadata_jp.json"))).T

    # already have start & end times (from before): start_time_sec & end_time_sec
    # now, will calculate real seizure times & adjusted times to 1970s

    real_start = pd.to_datetime(metadata.loc[subject_name]['Start day/time']).tz_convert(tz)
    #time delta between real start and 1970-01-01 00:00
    delta = pd.to_datetime(metadata.loc[subject_name]['Start day/time']).tz_convert(tz)-pd.to_datetime('1970-1-1').tz_localize(tz)
    delta_time  = delta-pd.to_timedelta(delta.days, unit='days')

    #start time with ref to real date but correct time of day
    sz_start_real= real_start + pd.to_timedelta(start_time_sec + 60, unit='s')  # real seizure start time does NOT go back 60 sec
    sz_end_real = real_start + pd.to_timedelta(end_time_sec, unit='s')

    #start time with ref to real date but correct time of day
    seg_start_real = real_start + pd.to_timedelta(start_time_sec, unit='s')
    seg_end_real = real_start + pd.to_timedelta(end_time_sec, unit='s')

    #start time with ref to 1970 date but correct time of day
    seg_start_1970= seg_start_real - pd.to_timedelta(delta.days, unit='days')
    seg_end_1970 = seg_end_real- pd.to_timedelta(delta.days, unit='days')


    #get data from portal
    portal_name = metadata.loc[metadata_subjectname]['portal_ID']
    ds = s.open_dataset(portal_name)


    ### Get EEG Data:

    data_pull_min=1
    clip_duration_sec =  end_time_sec - start_time_sec
    clip_duration_min = clip_duration_sec / 60

    # how many data_pull_min minute data pulls are there?
    n_iter = int(np.ceil(clip_duration_min / data_pull_min))


    eeg_channels = ['LA01',
    'LA02',
    'LA03']

    ### Read Saved Data:

    iEEG_filename = iEEG_filename
    file = f"../data/eeg_pt-{iEEG_filename}_start-{start_time_sec}_end-{end_time_sec}.h5" %()
    eeg_all = pd.read_hdf(file, key='ieeg')

    # shift to real time
    eeg_all.index = eeg_all.index + pd.to_datetime('1970-1-1').tz_localize(tz) + delta

    # define the global varibale so that it can be called in the plot function afterwards
 #   global eeg_channels

    eeg_channels = eeg_all.columns[[channel_num1 - 1, channel_num2 - 1, channel_num3 - 1, channel_num4 - 1]]

    eeg_channel_signals = signal_ref[[eeg_channels[0], eeg_channels[1], eeg_channels[2], eeg_channels[3]]]


    ### Shifting Time & Cut Data:

    #fetch saved ecg hr, watch and acc data
    data_root = "/gdrive/public/DATA/Human_Data/LB3_PIONEER/"
    ecg_hr = fetch_h5(subject_name, data_root, "ecg_heart")
    watch_hr = fetch_h5(subject_name, data_root, "watch_heart")
    watch_acc = fetch_h5(subject_name, data_root, "watch_acc")


    # #shift ecg/wearables data to 1970 keeping time correct
    # ecg_hr.index= ecg_hr.index - pd.to_timedelta(delta.days, unit='days')
    # watch_hr.index= watch_hr.index - pd.to_timedelta(delta.days, unit='days')
    # watch_acc.index= watch_acc.index - pd.to_timedelta(delta.days, unit='days')

    # #cutting data to selected range
    ecg_hr = ecg_hr[(ecg_hr.index>seg_start_real) & (ecg_hr.index< seg_end_real)]
    watch_hr = watch_hr[(watch_hr.index>seg_start_real) & (watch_hr.index< seg_end_real)]
    watch_acc = watch_acc[(watch_acc.index>seg_start_real) & (watch_acc.index< seg_end_real)]

    # Calc Acc magnitude
    watch_acc['mag'] = np.sqrt(np.sum(watch_acc**2, axis=1))

    ### LINE BELOW IS FOR DEBUGGING: ----- returning size & length of ECG, watch_hr, watch_acc
    # return sys.getsizeof(ecg_hr), sys.getsizeof(watch_hr), sys.getsizeof(watch_acc), len(ecg_hr), len(watch_hr), len(watch_acc)

    #########################################

    ### ### THIRD, creating 2 new separate data frames, one for the ECG signal & one for the watch/wearable data:

    # already have the EEG signal ready from before as the variable: signal_ref (must index specific channels of interest)

    # format resulting data into pandas DataFrame -- ECG
    ecg_data = {'ECG_HR': ecg_hr['heartRate']}
    ecg_df = pd.DataFrame(ecg_data)
    # ecg_df.index = pd.to_timedelta(t_sec, unit="S")

    # format resulting data into pandas DataFrame -- watch/wearable data
    watch_data = {'Watch_HR': watch_hr['heartRate'], 'Watch_Acc': watch_acc['mag']}
    watch_df = pd.DataFrame(watch_data)
    # watch_df.index = pd.to_timedelta(t_sec, unit="S")
 

    # returns the 4 specific channels of interest from ther EEG (instead of ALL the channels)...
    # FOLLOWING LINE IS FOR DEBUGGING:
    # return watch_hr
    return 'EEG Signal:', eeg_channel_signals, 'ECG Signal:', ecg_df, 'Watch HR and Accelerometer:', watch_df

    #########################################

    ### ### TRYING TO CONCATENATE AND JOIN THE EEG SIGNAL, ECG SIGNAL, 7 ALL WATCH WEARABLE DATA ONTO A SINGLE DATA FRAME:

    ### METHOD 1:
    # signal_filt vs. signal_ref        [signal_ref here)]
    # eeg_wearabale_dataFrame = signal_ref.assign(ECG_HR=ecg_hr, Watch_HR=watch_hr, Watch_Acc=watch_acc['mag'])

    ### METHOD 2:
    # signal_filt vs. signal_ref        [signal_ref here]
    # eeg_wearabale_dataFrame = signal_ref   
    ## ECG_HR data column added to data frame
    # eeg_wearabale_dataFrame['ECG_HR'] = ecg_hr
    ## Watch_HR data column added to data frame
    # eeg_wearabale_dataFrame['Watch_HR'] = watch_hr
    ## Watch_Acc data column added to data frame
    # eeg_wearabale_dataFrame['Watch_Acc'] = watch_acc['mag']

    ### METHOD 3:
    # signal_filt vs. signal_ref        [signal_ref here]
    # eeg_wearabale_dataFrame = signal_ref
    # eeg_wearabale_dataFrame.insert(-1, "ECG_HR", ecg_hr, True)
    # eeg_wearabale_dataFrame.insert(-1, "Watch_HR", watch_hr, True)
    # eeg_wearabale_dataFrame.insert(-1, "Watch_Acc", watch_acc['mag'], True)


    # # # CODE FOR DEBUGGING # # # 
 #   print(eeg_wearabale_dataFrame)
 #   if eeg_wearabale_dataFrame is None:
 #       return 'None'
 #   else: 
 #       return eeg_wearabale_dataFrame
    # # #                    # # # 

 #   eeg_wearabale_dataFrame.index = pd.to_timedelta(t_sec, unit="S")


    # create an hdf file out of the combined EEG & wearable data:
    # (from earlier): signal_ref.to_hdf(f"../data/eeg_pt-{iEEG_filename}_start-{start_time_sec}_end-{end_time_sec}.h5", key='ieeg')  

 #   return eeg_wearable_ref

    ### ### 
    ### ### OR, create a new DataFrame including ONLY the selected 4 EEG channels (eeg_channels variable)
    ### ### 

    # df = pd.DataFrame(EEG_ch1=eeg_channels[0])
    # eeg4_wearabale_dataFrame = df.assign(EEG_ch2=eeg_channels[1], EEG_ch3=eeg_channels[2], EEG_ch4=eeg_channels[3], ECG_HR=ecg_hr, Watch_HR=watch_hr, Watch_Acc=watch_acc)
    # eeg4_wearabale_dataFrame.index = pd.to_timedelta(t_sec, unit="S")

    # eeg4_wearable_ref = _bipolar_reference(eeg4_wearabale_dataFrame)

    # create an hdf file out of the combined EEG & wearable data:
    # (from earlier): signal_ref.to_hdf(f"../data/eeg_pt-{iEEG_filename}_start-{start_time_sec}_end-{end_time_sec}.h5", key='ieeg') 

    #  return eeg4_wearable_ref




