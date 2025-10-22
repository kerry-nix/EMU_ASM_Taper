
import sys
import os
root_folder = "./../"

import getopt
import numpy as np
import multiprocessing # todo
from os.path import join as ospj
import sys
import pandas as pd
from ieeg.auth import Session
import matplotlib.pyplot as plt
import matplotlib
from tqdm import tqdm
from tools import get_iEEG_data, bandpower
from numpy.lib.stride_tricks import sliding_window_view
from scipy.signal import iirnotch, filtfilt, butter, lfilter, medfilt
tz = 'US/Eastern'
import seaborn as sns

import warnings
warnings.filterwarnings("ignore")

import biosppy
from biosppy.signals import ecg
import pyhrv
import pyhrv.tools as tools
from pyhrv.hrv import hrv
import neurokit2 as nk
import heartpy as hp
import pyhrv.time_domain as td 
import pyhrv.frequency_domain as fd

from datetime import datetime, date, timedelta
import importlib


from scipy import signal
from scipy.signal import butter, lfilter, freqz, filtfilt, find_peaks, medfilt, ellip
from scipy.signal import savgol_filter, hilbert
import json


username = 'jalpanchal'
pwd_bin_path = "/gdrive/public/USERS/jalp/tools/jal_ieeglogin.bin"
with open(pwd_bin_path, "r") as f:
    s = Session(username, f.read())   



def offsetTo1970(data_raw, delta):
    tz = 'US/Eastern'
    data = data_raw.copy()
    data.index = data.index.tz_convert(tz)-timedelta(days=delta.days)
    return data


def fetch_h5(subject_name, root_folder, data_select):
    #root folder contains the LB3_oox_phasey subject folder
    #subject name : LB3_XXX_phaseY
    h5_files = [c for c in os.listdir(root_folder+subject_name+"/"+"wearables/pre-processed/") if (c.startswith(data_select)) & c.endswith('.h5') ]
    print(h5_files)
    h5_data = pd.read_hdf(root_folder+subject_name+"/"+"wearables/pre-processed/"+h5_files[0])
    return h5_data


def roliing_windows(start_s, end_s, winsize_s, overlap_s, fs):
    win_idx = []
    start_idx = start_s*fs
    end_idx = end_s*fs
    winsize_idx = winsize_s*fs
    overlap_idx = overlap_s*fs

    start_idx_all = np.arange(start_idx, end_idx-winsize_idx,winsize_idx-overlap_idx).reshape([-1,1])
    end_idx_all = np.arange(start_idx + winsize_idx, end_idx, winsize_idx-overlap_idx).reshape([-1,1])
    return np.concatenate([start_idx_all, end_idx_all], axis=1)

def calc_windows(start, stop, win_width_s, overlap_s):
    windows = []
    start = pd.to_datetime(start, unit='s')
    stop = pd.to_datetime(stop, unit='s')
    while start + pd.to_timedelta(win_width_s, unit='s') < stop:
        windows.append([start, start + pd.to_timedelta(win_width_s, unit='s')])
        start += pd.to_timedelta(win_width_s-overlap_s, unit='s')

    windows.append([stop -pd.to_timedelta(win_width_s, unit='s'), stop])    
    
    return windows


def plot_summary(hr_data_mean, hr_data_std, acc_data_mean, acc_data_std, ecg_hr_mean, ecg_hr_std, subject_name, str_time, end_time, sz_start=None, sleep_data=None):    

    #plot HR for all days
    plt.rcParams.update({'font.size': 16})

    
    fig,ax = plt.subplots(6,1, figsize=(15,25), sharex=True)    
    
    #HR plot
    axi = plt.subplot(6,1,1)
    axi.plot(hr_data_mean.index, hr_data_mean.heartRate, color='C3', alpha = 1)
    axi.fill_between(hr_data_mean.index, hr_data_mean.heartRate-hr_data_std.heartRate, hr_data_mean.heartRate + hr_data_std.heartRate, alpha = 0.3, facecolor='C3')            
    axi.set_ylabel("HR watch [bpm]")
    axi.grid(alpha=0.4)
    axi.set_title("Subject %s" %(subject_name))
    for l in sz_start:
            axi.axvline(l, color='C3')


    #HR plot ECG
    axi = plt.subplot(6,1,2)
    axi.plot(ecg_hr_mean.index, ecg_hr_mean.heartRate, color='C3', alpha = 1)
    axi.fill_between(ecg_hr_mean.index, ecg_hr_mean.heartRate-ecg_hr_std.heartRate, ecg_hr_mean.heartRate + ecg_hr_std.heartRate, alpha = 0.3, facecolor='C3')            
    axi.set_ylabel("HR ECG [bpm]")
    axi.grid(alpha=0.4)
    for l in sz_start:
            axi.axvline(l, color='C3')

    #rmssd
    axi = plt.subplot(6,1,3)
    axi.plot(ecg_hr_mean.index, ecg_hr_mean.rmssd, color='C4', alpha = 1)
    axi.fill_between(ecg_hr_mean.index, ecg_hr_mean.rmssd-ecg_hr_std.rmssd, ecg_hr_mean.rmssd + ecg_hr_std.rmssd, alpha = 0.3, facecolor='C4')            
    axi.set_ylabel("RMSSD [ms]")
    axi.grid(alpha=0.4)
    for l in sz_start:
            axi.axvline(l, color='C3')

    #hf
    axi = plt.subplot(6,1,4)
    axi.plot(ecg_hr_mean.index, ecg_hr_mean.hf, color='C4', alpha = 1)
    axi.fill_between(ecg_hr_mean.index, ecg_hr_mean.hf-ecg_hr_std.hf, ecg_hr_mean.hf + ecg_hr_std.hf, alpha = 0.3, facecolor='C4')            
    axi.set_ylabel("hf [%]")
    axi.grid(alpha=0.4)
    for l in sz_start:
            axi.axvline(l, color='C3')
    
    #Acc Plot

    #x axis
    axi = plt.subplot(6,1,5)
    axi.plot(acc_data_mean.index, acc_data_mean, color='C0', alpha = 1)
    axi.fill_between(acc_data_mean.index, acc_data_mean-acc_data_std, acc_data_mean + acc_data_std, alpha = 0.3, facecolor='C1')            
    axi.set_ylabel("Acc Mag[G]")
    axi.grid(alpha=0.4)
   

    # for a in ax:
    #     if sz_start.any():
    #         for l in sz_start:
    #             a.axvline(l, color='C3')
    for axi in ax:
        ylim = axi.get_ylim()
        if np.any(sleep_data):
            for sl in sleep_data.index:
                a.fill_betweenx(ylim, sleep_data.loc[sl,'startDate'], sleep_data.loc[sl,'endDate'], alpha=0.2, color ='black')



    plt.xlim([str_time, end_time])
    plt.gcf().autofmt_xdate() 
    
    plt.show()