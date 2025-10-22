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

### ### ###
# the following class is based off of & uses the 2 functions made under tools:
        # eeg_wearable.py
        # plot_eeg_wearable.py
### ### ###


class eeg_wearable_df_plot:

    def __init__(self, iEEG_filename, start_time, duration, subject_name, metadata_subjectname, channel_num1, channel_num2, channel_num3, channel_num4, yn_ecg, yn_watchHR, yn_watchAcc):

        self.iEEG_filename = iEEG_filename
        # time in seconds...
        self.start_time = start_time
        self.duration = duration
        self.subject_name = subject_name
        # subject_name & metadata_subjectname are often the same ; only differnet when the recording was interrupted, resulting in several snippets of the recording 
        self.metadata_subjectname = metadata_subjectname

        # channel number = NOT in Python indexing (starting from 1...)
        self.channel_num1 = channel_num1
        self.channel_num2 = channel_num2
        self.channel_num3 = channel_num3
        self.channel_num4 = channel_num4

        # for the following three:
            # 1 = the data file exists in the pre-processed data folder for the subject
            # 0 = the data file does NOT exist in the pre-processed data folder for the subject
        self.yn_ecg = yn_ecg
        self.yn_watchHR = yn_watchHR
        self.yn_watchAcc = yn_watchAcc



        # variables to eventually be declared as properties/objects:
            # ... eeg_channel_signals
            # ... sz_start_real
            # ... ecg_df
            # ... watch_df
            # ... eeg_channels
            # ... signal_ref
            # ... eeg_all

    # function for DISPLAYING the dataframes + EEG signals

    def display_eeg_wearable(self):

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
        start_time_sec = self.start_time - 60    # showing starting 1 minute before the seizure
        end_time_sec = self.start_time + self.duration

        # accessing iEEG file name & specific time clip:
        iEEG_filename = self.iEEG_filename
        start_time_usec = start_time_sec * 1e6
        end_time_usec = end_time_sec * 1e6

        # obtaining iEEG data
        data, fs = get_iEEG_data(username, pwd_bin_file, self.iEEG_filename, start_time_usec, end_time_usec, select_electrodes= None)

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

        self.signal_ref = _bipolar_reference(signal_filt)

        ### return signal_ref ###

        # create an hdf file out of the EEG data 
        #signal_ref.to_hdf(f"../data/eeg_pt-{iEEG_filename}_start-{start_time_sec}_end-{end_time_sec}.h5", key='ieeg')

        # plotting the EEG channels (ONLY) to observe which channels show more seizure activity 
        plt.figure(figsize=(10, 100)); plt.plot(self.signal_ref.iloc[:, 0:self.signal_ref.shape[1]] + np.arange(self.signal_ref.shape[1])*-1000)

        #########################################

        ### ### SECOND: overall, to obtain the wearables data & EEG data from the step above to plot


        # obtaining subject name & data
        # subject_name = subject_name
        metadata = pd.DataFrame(json.load(open("./tools/subject_metadata_jp.json"))).T

        # already have start & end times (from before): start_time_sec & end_time_sec
        # now, will calculate real seizure times & adjusted times to 1970s

        real_start = pd.to_datetime(metadata.loc[self.subject_name]['Start day/time']).tz_convert(tz)
        #time delta between real start and 1970-01-01 00:00
        delta = pd.to_datetime(metadata.loc[self.subject_name]['Start day/time']).tz_convert(tz)-pd.to_datetime('1970-1-1').tz_localize(tz)
        delta_time  = delta-pd.to_timedelta(delta.days, unit='days')

        #start time with ref to real date but correct time of day
        self.sz_start_real= real_start + pd.to_timedelta(start_time_sec + 60, unit='s')  # real seizure start time does NOT go back 60 sec
        sz_end_real = real_start + pd.to_timedelta(end_time_sec, unit='s')

        #start time with ref to real date but correct time of day
        seg_start_real = real_start + pd.to_timedelta(start_time_sec, unit='s')
        seg_end_real = real_start + pd.to_timedelta(end_time_sec, unit='s')

        #start time with ref to 1970 date but correct time of day
        seg_start_1970= seg_start_real - pd.to_timedelta(delta.days, unit='days')
        seg_end_1970 = seg_end_real- pd.to_timedelta(delta.days, unit='days')


        #get data from portal
        portal_name = metadata.loc[self.metadata_subjectname]['portal_ID']
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

        # iEEG_filename = iEEG_filename
      #  file = f"../data/eeg_pt-{self.iEEG_filename}_start-{start_time_sec}_end-{end_time_sec}.h5" %()
      #  self.eeg_all = pd.read_hdf(file, key='ieeg')

        self.eeg_all = self.signal_ref

        # shift to real time
        self.eeg_all.index = self.eeg_all.index + pd.to_datetime('1970-1-1').tz_localize(tz) + delta

        # define the global varibale so that it can be called in the plot function afterwards
     #   global eeg_channels

        self.eeg_channels = self.eeg_all.columns[[self.channel_num1 - 1, self.channel_num2 - 1, self.channel_num3 - 1, self.channel_num4 - 1]]

        self.eeg_channel_signals = self.signal_ref[[self.eeg_channels[0], self.eeg_channels[1], self.eeg_channels[2], self.eeg_channels[3]]]
        # INDEX = ??? --> using index of eeg_all (which should be the same as the index for this)


        ### Shifting Time & Cut Data:

        #fetch saved ecg hr, watch and acc data
        data_root = "/gdrive/public/DATA/Human_Data/LB3_PIONEER/"

        # ONLY to fetch the data IF the data file exists (some subjects do NOT have certain data)
        if self.yn_ecg == 1:
            ecg_hr = fetch_h5(self.subject_name, data_root, "ecg_heart")
        if self.yn_watchHR == 1:
            watch_hr = fetch_h5(self.subject_name, data_root, "watch_heart")
        if self.yn_watchAcc == 1:
            watch_acc = fetch_h5(self.subject_name, data_root, "watch_acc")


        # #shift ecg/wearables data to 1970 keeping time correct
        # ecg_hr.index= ecg_hr.index - pd.to_timedelta(delta.days, unit='days')
        # watch_hr.index= watch_hr.index - pd.to_timedelta(delta.days, unit='days')
        # watch_acc.index= watch_acc.index - pd.to_timedelta(delta.days, unit='days')


        # # cutting data to selected range + creating data frame:
                # data only added if the data file exists (1 = exists ; 0 = does NOT exist)
    
        if self.yn_ecg == 1:
            ecg_hr = ecg_hr[(ecg_hr.index>seg_start_real) & (ecg_hr.index< seg_end_real)]
            # format resulting data into pandas DataFrame -- ECG
            ecg_data = {'ECG_HR': ecg_hr['heartRate']}
            self.ecg_df = pd.DataFrame(ecg_data)
         #   self.ecg_df.index = pd.to_timedelta(t_sec, unit="S")



        # print('TEST')
        if self.yn_watchHR == 1 and self.yn_watchAcc == 1:
            watch_hr = watch_hr[(watch_hr.index>seg_start_real) & (watch_hr.index< seg_end_real)]
            watch_acc = watch_acc[(watch_acc.index>seg_start_real) & (watch_acc.index< seg_end_real)]
            watch_acc['mag'] = np.sqrt(np.sum(watch_acc**2, axis=1))
            # format resulting data into pandas DataFrame -- watch data
            watchHR_data = {'Watch_HR': watch_hr['heartRate']}
            watchAcc_data = {'Watch_Acc': watch_acc['mag']}
            # watch_data = {'Watch_HR': watch_hr['heartRate'], 'Watch_Acc': watch_acc['mag']}
            #
            self.watchHR_df = pd.DataFrame(watchHR_data)
            self.watchAcc_df = pd.DataFrame(watchAcc_data)
            # print(self.watch_df.index.is_unique) # --> True
            # self.watch_df['Watch_Acc'] = watch_acc['mag']

        elif self.yn_watchHR == 1 and self.yn_watchAcc != 1:
            watch_hr = watch_hr[(watch_hr.index>seg_start_real) & (watch_hr.index< seg_end_real)]
            # dict --> df
            watchHR_data = {'Watch_HR': watch_hr['heartRate']}
            ###
            # temp = {'Watch_Acc': 'No Watch Accelerometer Data'}
            ###
            temp = {'Watch_Acc': [0]*len(watch_hr['heartRate'])}
            # format resulting data into pandas DataFrame -- watch data
            self.watchHR_df = pd.DataFrame(watchHR_data)
            self.watchAcc_df = pd.DataFrame(temp)
            
        elif self.yn_watchHR != 1 and self.yn_watchAcc == 1:
            watch_acc = watch_acc[(watch_acc.index>seg_start_real) & (watch_acc.index< seg_end_real)]
            watch_acc['mag'] = np.sqrt(np.sum(watch_acc**2, axis=1))
            # dict --> df
            watchAcc_data = {'Watch_Acc': watch_acc['mag']}
            ###
            # temp = {'Watch_HR': 'No Watch HR Data'}
            ###
            temp = {'Watch_HR': [0]*len(watch_acc['mag'])}
            # format resulting data into pandas DataFrame -- watch data
            self.watchAcc_df = pd.DataFrame(watchAcc_data)
            self.watchHR_df = pd.DataFrame(temp)
        # else:
        #     watch_data = 0
        
        # # converting into data frame (for watch):
        # if watch_data != 0:
        #     self.watch_df = pd.DataFrame(watch_data)
           # #   self.watch_df.index = pd.to_timedelta(t_sec, unit="S")
        
        #################################################################################################

        # if self.yn_watchHR == 1:
        #     watch_hr = watch_hr[(watch_hr.index>seg_start_real) & (watch_hr.index< seg_end_real)]
        #     # format resulting data into pandas DataFrame -- watch_HR
        #     watch_data = {'Watch_HR': watch_hr['heartRate']}
        # if self.yn_watchAcc == 1:
        #     watch_acc = watch_acc[(watch_acc.index>seg_start_real) & (watch_acc.index< seg_end_real)]
        #     # Calc Acc magnitude
        #     watch_acc['mag'] = np.sqrt(np.sum(watch_acc**2, axis=1))
        #     # format resulting data into pandas DataFrame -- watch_Acc
        #     if self.yn_watchHR == 1:
        #         watch_data['Watch_Acc'] = watch_acc['mag']      # add a column
        #     else:
        #         watch_data = {'Watch_Acc': watch_acc['mag']}    # create the dataframe
       
        # # converting into data frame (for watch):
        # if self.yn_watchHR or self.yn_watchAcc == 1:
        #     self.watch_df = pd.DataFrame(watch_data)
        #     #   self.watch_df.index = pd.to_timedelta(t_sec, unit="S")

        #################################################################################################


        ### LINE BELOW IS FOR DEBUGGING: ----- returning size & length of ECG, watch_hr, watch_acc
        # return sys.getsizeof(ecg_hr), sys.getsizeof(watch_hr), sys.getsizeof(watch_acc), len(ecg_hr), len(watch_hr), len(watch_acc)

        #########################################

        ### ### THIRD, creating 2 new separate data frames, one for the ECG signal & one for the watch/wearable data:

        # already have the EEG signal ready from before as the variable: signal_ref (must index specific channels of interest)

    

        # returns the 4 specific channels of interest from ther EEG (instead of ALL the channels)...
        # FOLLOWING LINE IS FOR DEBUGGING:
        # return watch_hr

        ### DEBUGGING CODE ###
     #   return len(self.ecg_df), self.ecg_df.size, len(self.watch_df), self.watch_df.size
     #   return len(ecg_hr.index), len(ecg_hr['heartRate'])
     #   return len(self.ecg_df.index), len(self.watch_df.index)
        ###                ### 

        return 'EEG Signal:', self.eeg_channel_signals, 'ECG Signal:', self.ecg_df, 'Watch HR:', self.watchHR_df, 'Watch Accelerometer:', self.watchAcc_df


 ############# 

 # function for SAVING the dataframes

    def save_eeg_wearable(self):
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
            start_time_sec = self.start_time - 60    # showing starting 1 minute before the seizure
            end_time_sec = self.start_time + self.duration

            # accessing iEEG file name & specific time clip:
            iEEG_filename = self.iEEG_filename
            start_time_usec = start_time_sec * 1e6
            end_time_usec = end_time_sec * 1e6

            # obtaining iEEG data
            data, fs = get_iEEG_data(username, pwd_bin_file, self.iEEG_filename, start_time_usec, end_time_usec, select_electrodes= None)

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

            self.signal_ref = _bipolar_reference(signal_filt)

            ### return signal_ref ###

            # create an hdf file out of the EEG data 
            #signal_ref.to_hdf(f"../data/eeg_pt-{iEEG_filename}_start-{start_time_sec}_end-{end_time_sec}.h5", key='ieeg')

            # plotting the EEG channels (ONLY) to observe which channels show more seizure activity 
            # plt.figure(figsize=(10, 100)); plt.plot(self.signal_ref.iloc[:, 0:self.signal_ref.shape[1]] + np.arange(self.signal_ref.shape[1])*-1000)

            #########################################

            ### ### SECOND: overall, to obtain the wearables data & EEG data from the step above to plot


            # obtaining subject name & data
            # subject_name = subject_name
            metadata = pd.DataFrame(json.load(open("./tools/subject_metadata_jp.json"))).T

            # already have start & end times (from before): start_time_sec & end_time_sec
            # now, will calculate real seizure times & adjusted times to 1970s

            real_start = pd.to_datetime(metadata.loc[self.subject_name]['Start day/time']).tz_convert(tz)
            #time delta between real start and 1970-01-01 00:00
            delta = pd.to_datetime(metadata.loc[self.subject_name]['Start day/time']).tz_convert(tz)-pd.to_datetime('1970-1-1').tz_localize(tz)
            delta_time  = delta-pd.to_timedelta(delta.days, unit='days')

            #start time with ref to real date but correct time of day
            self.sz_start_real= real_start + pd.to_timedelta(start_time_sec + 60, unit='s')  # real seizure start time does NOT go back 60 sec
            sz_end_real = real_start + pd.to_timedelta(end_time_sec, unit='s')

            #start time with ref to real date but correct time of day
            seg_start_real = real_start + pd.to_timedelta(start_time_sec, unit='s')
            seg_end_real = real_start + pd.to_timedelta(end_time_sec, unit='s')

            #start time with ref to 1970 date but correct time of day
            seg_start_1970= seg_start_real - pd.to_timedelta(delta.days, unit='days')
            seg_end_1970 = seg_end_real- pd.to_timedelta(delta.days, unit='days')


            #get data from portal
            portal_name = metadata.loc[self.metadata_subjectname]['portal_ID']
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

            # iEEG_filename = iEEG_filename
        #  file = f"../data/eeg_pt-{self.iEEG_filename}_start-{start_time_sec}_end-{end_time_sec}.h5" %()
        #  self.eeg_all = pd.read_hdf(file, key='ieeg')

            self.eeg_all = self.signal_ref

            # shift to real time
            self.eeg_all.index = self.eeg_all.index + pd.to_datetime('1970-1-1').tz_localize(tz) + delta

            # define the global varibale so that it can be called in the plot function afterwards
        #   global eeg_channels

            self.eeg_channels = self.eeg_all.columns[[self.channel_num1 - 1, self.channel_num2 - 1, self.channel_num3 - 1, self.channel_num4 - 1]]

            self.eeg_channel_signals = self.signal_ref[[self.eeg_channels[0], self.eeg_channels[1], self.eeg_channels[2], self.eeg_channels[3]]]
            # INDEX = ??? --> using index of eeg_all (which should be the same as the index for this)


            ### Shifting Time & Cut Data:

            #fetch saved ecg hr, watch and acc data
            data_root = "/gdrive/public/DATA/Human_Data/LB3_PIONEER/"

            # ONLY to fetch the data IF the data file exists (some subjects do NOT have certain data)
            if self.yn_ecg == 1:
                ecg_hr = fetch_h5(self.subject_name, data_root, "ecg_heart")
            if self.yn_watchHR == 1:
                watch_hr = fetch_h5(self.subject_name, data_root, "watch_heart")
            if self.yn_watchAcc == 1:
                watch_acc = fetch_h5(self.subject_name, data_root, "watch_acc")


            # #shift ecg/wearables data to 1970 keeping time correct
            # ecg_hr.index= ecg_hr.index - pd.to_timedelta(delta.days, unit='days')
            # watch_hr.index= watch_hr.index - pd.to_timedelta(delta.days, unit='days')
            # watch_acc.index= watch_acc.index - pd.to_timedelta(delta.days, unit='days')


            # # cutting data to selected range + creating data frame:
                    # data only added if the data file exists (1 = exists ; 0 = does NOT exist)
        
            if self.yn_ecg == 1:
                ecg_hr = ecg_hr[(ecg_hr.index>seg_start_real) & (ecg_hr.index< seg_end_real)]
                # format resulting data into pandas DataFrame -- ECG
                ecg_data = {'ECG_HR': ecg_hr['heartRate']}
                self.ecg_df = pd.DataFrame(ecg_data)
            #   self.ecg_df.index = pd.to_timedelta(t_sec, unit="S")



            # print('TEST')
            if self.yn_watchHR == 1 and self.yn_watchAcc == 1:
                watch_hr = watch_hr[(watch_hr.index>seg_start_real) & (watch_hr.index< seg_end_real)]
                watch_acc = watch_acc[(watch_acc.index>seg_start_real) & (watch_acc.index< seg_end_real)]
                watch_acc['mag'] = np.sqrt(np.sum(watch_acc**2, axis=1))
                # format resulting data into pandas DataFrame -- watch data
                watchHR_data = {'Watch_HR': watch_hr['heartRate']}
                watchAcc_data = {'Watch_Acc': watch_acc['mag']}
                # watch_data = {'Watch_HR': watch_hr['heartRate'], 'Watch_Acc': watch_acc['mag']}
                #
                self.watchHR_df = pd.DataFrame(watchHR_data)
                self.watchAcc_df = pd.DataFrame(watchAcc_data)
                # print(self.watch_df.index.is_unique) # --> True
                # self.watch_df['Watch_Acc'] = watch_acc['mag']

            elif self.yn_watchHR == 1 and self.yn_watchAcc != 1:
                watch_hr = watch_hr[(watch_hr.index>seg_start_real) & (watch_hr.index< seg_end_real)]
                # dict --> df
                watchHR_data = {'Watch_HR': watch_hr['heartRate']}
                ###
                # temp = {'Watch_Acc': 'No Watch Accelerometer Data'}
                ###
                temp = {'Watch_Acc': [0]*len(watch_hr['heartRate'])}
                # format resulting data into pandas DataFrame -- watch data
                self.watchHR_df = pd.DataFrame(watchHR_data)
                self.watchAcc_df = pd.DataFrame(temp)
                
            elif self.yn_watchHR != 1 and self.yn_watchAcc == 1:
                watch_acc = watch_acc[(watch_acc.index>seg_start_real) & (watch_acc.index< seg_end_real)]
                watch_acc['mag'] = np.sqrt(np.sum(watch_acc**2, axis=1))
                # dict --> df
                watchAcc_data = {'Watch_Acc': watch_acc['mag']}
                ###
                # temp = {'Watch_HR': 'No Watch HR Data'}
                ###
                temp = {'Watch_HR': [0]*len(watch_acc['mag'])}
                # format resulting data into pandas DataFrame -- watch data
                self.watchAcc_df = pd.DataFrame(watchAcc_data)
                self.watchHR_df = pd.DataFrame(temp)
            # else:
            #     watch_data = 0
            
            # # converting into data frame (for watch):
            # if watch_data != 0:
            #     self.watch_df = pd.DataFrame(watch_data)
            # #   self.watch_df.index = pd.to_timedelta(t_sec, unit="S")
            
            #################################################################################################

            # if self.yn_watchHR == 1:
            #     watch_hr = watch_hr[(watch_hr.index>seg_start_real) & (watch_hr.index< seg_end_real)]
            #     # format resulting data into pandas DataFrame -- watch_HR
            #     watch_data = {'Watch_HR': watch_hr['heartRate']}
            # if self.yn_watchAcc == 1:
            #     watch_acc = watch_acc[(watch_acc.index>seg_start_real) & (watch_acc.index< seg_end_real)]
            #     # Calc Acc magnitude
            #     watch_acc['mag'] = np.sqrt(np.sum(watch_acc**2, axis=1))
            #     # format resulting data into pandas DataFrame -- watch_Acc
            #     if self.yn_watchHR == 1:
            #         watch_data['Watch_Acc'] = watch_acc['mag']      # add a column
            #     else:
            #         watch_data = {'Watch_Acc': watch_acc['mag']}    # create the dataframe
        
            # # converting into data frame (for watch):
            # if self.yn_watchHR or self.yn_watchAcc == 1:
            #     self.watch_df = pd.DataFrame(watch_data)
            #     #   self.watch_df.index = pd.to_timedelta(t_sec, unit="S")

            #################################################################################################


            ### LINE BELOW IS FOR DEBUGGING: ----- returning size & length of ECG, watch_hr, watch_acc
            # return sys.getsizeof(ecg_hr), sys.getsizeof(watch_hr), sys.getsizeof(watch_acc), len(ecg_hr), len(watch_hr), len(watch_acc)

            #########################################

            ### ### THIRD, creating 2 new separate data frames, one for the ECG signal & one for the watch/wearable data:

            # already have the EEG signal ready from before as the variable: signal_ref (must index specific channels of interest)

        

            # returns the 4 specific channels of interest from ther EEG (instead of ALL the channels)...
            # FOLLOWING LINE IS FOR DEBUGGING:
            # return watch_hr

            ### DEBUGGING CODE ###
        #   return len(self.ecg_df), self.ecg_df.size, len(self.watch_df), self.watch_df.size
        #   return len(ecg_hr.index), len(ecg_hr['heartRate'])
        #   return len(self.ecg_df.index), len(self.watch_df.index)
            ###                ### 

            return self.eeg_channel_signals, self.ecg_df, self.watchHR_df, self.watchAcc_df

        ### NOW, function to plot...

    def plot_eeg_wearable(self):
   
        # using the 'eeg_channels", "ecg_hr", "watch_hr", "watch_acc", and "sz_start_real" variables from eeg_wearable.py function 

        # yn_ecg, yn_watchHR, yn_watchAcc

        #########################################

        ### ### FOURTH, plot the EEG & werable data:

        rw = len(self.eeg_channels)+3

        fig, ax = plt.subplots(rw, 1, figsize=(10,rw), sharex=True)
        sns.set_style('whitegrid')
        plt.rcParams.update({'font.size': 10})

        for i in range(len(self.eeg_channels)):
            axi = ax[i]
            axi.plot(self.eeg_all.index, self.eeg_channel_signals[self.eeg_channels[i]])
            axi.set_ylabel("%s\n[uV]" %(self.eeg_channels[i]))
        
        if self.yn_ecg == 1:
            axi = ax[len(self.eeg_channels)]
            axi.plot(self.ecg_df.index, self.ecg_df['ECG_HR'],'o-', color='C3')
            axi.set_ylabel("ECG HR\n[bpm]")
            axi.set_ylim([self.ecg_df['ECG_HR'].min() - 5, self.ecg_df['ECG_HR'].max() + 5])

        if self.yn_watchHR == 1:
            # counting to make sure the watch HR data does not contain NaN values
            countNaN_watchHR = self.watchHR_df['Watch_HR'].isna().sum()
            # # #
            if countNaN_watchHR == 0:
                if self.watchHR_df.empty:
                    pass
                else:
                    axi = ax[len(self.eeg_channels)+1]
                    axi.plot(self.watchHR_df['Watch_HR'].dropna(),'o-', color='C3')
                    axi.set_ylabel("Watch HR\n[bpm]")
                    axi.set_ylim([self.watchHR_df['Watch_HR'].min() - 5, self.watchHR_df['Watch_HR'].max() + 5])


        if self.yn_watchAcc == 1:
            # counting to make sure the watch Acc data does not contain NaN values
            countNaN_watchAcc = self.watchAcc_df['Watch_Acc'].isna().sum()
            # # #
            if countNaN_watchAcc == 0:
                if self.watchAcc_df.empty:
                    pass
                else:
                    axi = ax[len(self.eeg_channels)+2]
                    axi.plot(self.watchAcc_df.index, self.watchAcc_df['Watch_Acc'], color='C1')
                    axi.set_ylabel("Acc Mag\n[G]")
                    axi.set_xlabel("Time of day")
                    axi.set_ylim([self.watchAcc_df['Watch_Acc'].min() - 1, self.watchAcc_df['Watch_Acc'].max() + 1])

        for a in ax:
            a.axvline(x=self.sz_start_real, ls='--', color='grey')


        # xformatter = matplotlib.dates.DateFormatter('%H:%M')
        # axi.xaxis.set_major_formatter(xformatter)

        plt.suptitle("Multilodal seizure recording")
        sns.despine(bottom=True)


        plt.tight_layout()
        # plt.savefig("./../../figures/%s_%s_%s.png" %(subject_name,start_time_sec,end_time_sec), dpi = 300, bbox_inches='tight')
        plt.show()

        # return fig
        
        # end function (should show plot already)

