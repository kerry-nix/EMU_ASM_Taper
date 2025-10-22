   

import tools
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
from scipy.signal import iirnotch, filtfilt, butter

import sys
import os
sys.path.insert(1, os.path.join("./tools/"))
from helpers import *

# reload library
importlib.reload(sys.modules['helpers'])
from helpers import *

   
def plot_eeg_wearable(eeg_channels, ecg_hr, watch_hr, watch_acc, sz_start_real):
   
    # using the 'eeg_channels", "ecg_hr", "watch_hr", "watch_acc", and "sz_start_real" variables from eeg_wearable.py function 

    #########################################

    ### ### FOURTH, plot the EEG & werable data:

    rw = len(eeg_channels)+3

    fig, ax = plt.subplots(rw, 1, figsize=(10,rw), sharex=True)
    sns.set_style('whitegrid')
    plt.rcParams.update({'font.size': 10})

    for i in range(len(eeg_channels)):
        axi = ax[i]
        axi.plot(eeg_all.index, eeg_all[eeg_channels[i]])
        axi.set_ylabel("%s\n[uV]" %(eeg_channels[i]))
    

    axi = ax[len(eeg_channels)]
    axi.plot(ecg_hr.index, ecg_hr['heartRate'],'o-', color='C3')
    axi.set_ylabel("ECG HR\n[bpm]")
    axi.set_ylim([35,140])

    axi = ax[len(eeg_channels)+1]
    axi.plot(watch_hr['heartRate'].dropna(),'o-', color='C3')
    axi.set_ylabel("Watch HR\n[bpm]")
    axi.set_ylim([40,100])

    axi = ax[len(eeg_channels)+2]
    axi.plot(watch_acc.index, watch_acc['mag'], color='C1')
    axi.set_ylabel("Acc Mag\n[G]")
    axi.set_xlabel("Time of day")

    for a in ax:
        a.axvline(x=sz_start_real, ls='--', color='grey')


    # xformatter = matplotlib.dates.DateFormatter('%H:%M')
    # axi.xaxis.set_major_formatter(xformatter)

    plt.suptitle("Multilodal seizure recording")
    sns.despine(bottom=True)


    plt.tight_layout()
    # plt.savefig("./../../figures/%s_%s_%s.png" %(subject_name,start_time_sec,end_time_sec), dpi = 300, bbox_inches='tight')
    plt.show()

    return fig

