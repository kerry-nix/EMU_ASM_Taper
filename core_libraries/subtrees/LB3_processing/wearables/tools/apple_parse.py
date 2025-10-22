

from datetime import date, time, datetime, timedelta
import pytz
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from dateutil import parser
import os
import seaborn as sns
tz = "US/Eastern"


# functions
def offsetTo1970(data_raw, delta):
    tz = 'US/Eastern'
    data = data_raw.copy()
    data.index = data.index.tz_convert(tz)-timedelta(days=delta.days)
    return data


def fetch_watch_hr_data_h5(subject_name, root_folder):
    h5_files = [c for c in os.listdir(root_folder+subject_name+"/"+"wearables/pre-processed/") if c.endswith('watch_heart.h5')]
    heart_h5 = pd.read_hdf(root_folder+subject_name+"/"+"wearables/pre-processed/"+h5_files[0])
    return heart_h5
    

def fetch_acc_data_h5(subject_name, root_folder):
    h5_files = [c for c in os.listdir(root_folder+subject_name+"/"+"wearables/pre-processed/") if c.endswith('acc.h5')]
    acc_h5 = pd.read_hdf(root_folder+subject_name+"/"+"wearables/pre-processed/"+h5_files[0])
    return acc_h5

def fetch_h5(subject_name, root_folder, data_select):
    #root folder contains the LB3_oox_phasey subject folder
    #subject name : LB3_XXX_phaseY
    h5_files = [c for c in os.listdir(root_folder+subject_name+"/"+"wearables/pre-processed/") if (c.startswith(data_select)) & c.endswith('.h5') ]
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
    


def plot_all_hr_data(hr_data, win_size, subject_name):
    #plot HR for all days
    plt.rcParams.update({'font.size': 18})
#     start_date = datetime(2021,10,19)
#     date = start_date.strftime("%Y-%m-%d")
#     num_days = 7
#     end_date = start_date + timedelta(days=num_days)

    sel_data = hr_data

    sel_data_avg = sel_data.resample(win_size).mean()
    sel_data_std = sel_data.resample(win_size).std()

    fig = plt.figure(figsize=(15,10))

    ax = fig.add_subplot(211)
    ax.plot(sel_data.index, sel_data.heartRate, color='C3', alpha = 0.2)
    ax.plot(sel_data_avg.index, sel_data_avg.heartRate, color='C3', alpha = 1)
    plt.gcf().autofmt_xdate()
    ax.set_ylabel("Heart Rate [bpm]")
    ax.set_xlabel("Date-Time")
    ax.grid(alpha=0.4)
    ax.set_title('Subject : %s' %(subject_name))

#     ax1 = fig.add_subplot(212, sharex=ax)
#     ax1.plot(sel_data.index, 1/sel_data.samplingInterval, alpha = 0.2, color = 'C1')
#     ax1.plot(sel_data_avg.index, 1/sel_data_avg.samplingInterval, alpha = 1, color = 'C1')
#     plt.gcf().autofmt_xdate()
#     ax1.set_ylabel("Sampling rate (Hz)")
#     ax1.set_xlabel("Date-Time")
#     ax1.set_ylim((0,1))
#     ax1.grid(alpha=0.4)

    plt.show()


def plot_all_acc_data(acc_data, win_size, subject_name):
    #plot HR for all days
    plt.rcParams.update({'font.size': 18})
#     start_date = datetime(2021,10,19)
#     date = start_date.strftime("%Y-%m-%d")
#     num_days = 7
#     end_date = start_date + timedelta(days=num_days)

    sel_data = acc_data

    sel_data_avg = sel_data.resample(win_size).mean()
    sel_data_std = sel_data.resample(win_size).std()

    fig = plt.figure(figsize=(15,10))

    ax = fig.add_subplot(311)
    
    ax.plot(sel_data_avg.index, sel_data_avg.accX, color='C1', alpha = 1)
    ax.fill_between(sel_data_std.index, sel_data_avg.accX-sel_data_std.accX, sel_data_avg.accX + sel_data_std.accX, alpha = 0.3, facecolor='C1')
    plt.gcf().autofmt_xdate()
    ax.set_ylabel("Acc X [G]")
    ax.set_xlabel("Date-Time")
    ax.grid(alpha=0.4)
    ax.set_title('Subject : %s' %(subject_name))

    ax1 = fig.add_subplot(312, sharex=ax)
    ax1.plot(sel_data_avg.index, sel_data_avg.accY, color='C0', alpha = 1)
    ax1.fill_between(sel_data_std.index, sel_data_avg.accY-sel_data_std.accY, sel_data_avg.accY + sel_data_std.accY, alpha = 0.3, facecolor='C0')
    plt.gcf().autofmt_xdate()
    ax1.set_ylabel("Acc Y [G]")
    ax1.grid(alpha=0.4)
    
    ax2 = fig.add_subplot(313, sharex=ax)
    ax2.plot(sel_data_avg.index, sel_data_avg.accZ, color='C2', alpha = 1)
    ax2.fill_between(sel_data_std.index, sel_data_avg.accZ-sel_data_std.accZ, sel_data_avg.accZ+sel_data_std.accZ, alpha = 0.3, facecolor='C2')
    plt.gcf().autofmt_xdate()
    ax2.set_ylabel("Acc Z [G]")
    ax2.grid(alpha=0.4)
    
    plt.show()
    
    
def plot_hr_acc(hr_data_avg, hr_data_std, acc_data_avg, acc_data_std, subject_name, str_time, end_time):
    
    #plot HR for all days
    plt.rcParams.update({'font.size': 18})
    
    fig = plt.figure(figsize=(15,10))
    
    #HR plot
    ax = fig.add_subplot(411)
    ax.plot(hr_data_avg.index, hr_data_avg.heartRate, color='C3', alpha = 1)
    ax.fill_between(hr_data_avg.index, hr_data_avg.heartRate-hr_data_std.heartRate, hr_data_avg.heartRate + hr_data_std.heartRate, alpha = 0.3, facecolor='C3')
    plt.gcf().autofmt_xdate()
    ax.set_ylabel("Heart Rate [bpm]")
    ax.set_xlabel("Date-Time")
    ax.set_xlim([str_time, end_time])
    ax.grid(alpha=0.4)
    ax.set_title('Subject : %s' %(subject_name))
    
    #Acc Plot
    ax1 = fig.add_subplot(412,sharex=ax)
    
    ax1.plot(acc_data_avg.index, acc_data_avg.accX, color='C1', alpha = 1)
    ax1.fill_between(acc_data_std.index, acc_data_avg.accX-acc_data_std.accX, acc_data_avg.accX + acc_data_std.accX, alpha = 0.3, facecolor='C1')
    plt.gcf().autofmt_xdate()
    ax1.set_ylabel("Acc X [G]")
    ax1.grid(alpha=0.4)

    ax2 = fig.add_subplot(413, sharex=ax)
    ax2.plot(acc_data_avg.index, acc_data_avg.accY, color='C0', alpha = 1)
    ax2.fill_between(acc_data_std.index, acc_data_avg.accY-acc_data_std.accY, acc_data_avg.accY + acc_data_std.accY, alpha = 0.3, facecolor='C0')
    plt.gcf().autofmt_xdate()
    ax2.set_ylabel("Acc Y [G]")
    ax2.grid(alpha=0.4)
    
    ax3 = fig.add_subplot(414, sharex=ax)
    ax3.plot(acc_data_avg.index, acc_data_avg.accZ, color='C2', alpha = 1)
    ax3.fill_between(acc_data_std.index, acc_data_avg.accZ-acc_data_std.accZ, acc_data_avg.accZ+acc_data_std.accZ, alpha = 0.3, facecolor='C2')
    plt.gcf().autofmt_xdate()
    ax3.set_ylabel("Acc Z [G]")
    ax3.grid(alpha=0.4)
    
    plt.show()


def plot_hr_acc_ecg(hr_data_avg, hr_data_std, acc_data_avg, acc_data_std, ecg_hr_mean, ecg_hr_std, subject_name, str_time, end_time):
    
    #plot HR for all days
    plt.rcParams.update({'font.size': 18})
    
    fig = plt.figure(figsize=(15,25))
    
    #HR plot
    ax = fig.add_subplot(811)
    ax.plot(hr_data_avg.index, hr_data_avg.heartRate, color='C3', alpha = 1)
    ax.fill_between(hr_data_avg.index, hr_data_avg.heartRate-hr_data_std.heartRate, hr_data_avg.heartRate + hr_data_std.heartRate, alpha = 0.3, facecolor='C3')
    plt.gcf().autofmt_xdate()
    ax.set_ylabel("HR watch [bpm]")
    ax.set_xlabel("Date-Time")
    ax.set_xlim([str_time, end_time])
    ax.grid(alpha=0.4)
    ax.set_title('Subject : %s' %(subject_name))


    #HR plot ECG
    ax = fig.add_subplot(812)
    ax.plot(ecg_hr_mean.index, ecg_hr_mean.heartRate, color='C3', alpha = 1)
    ax.fill_between(ecg_hr_mean.index, ecg_hr_mean.heartRate-ecg_hr_std.heartRate, ecg_hr_mean.heartRate + ecg_hr_std.heartRate, alpha = 0.3, facecolor='C3')
    plt.gcf().autofmt_xdate()
    ax.set_ylabel("HR ECG [bpm]")
    ax.set_xlabel("Date-Time")
    ax.set_xlim([str_time, end_time])
    ax.grid(alpha=0.4)

    #rmssd
    ax = fig.add_subplot(813)
    ax.plot(ecg_hr_mean.index, ecg_hr_mean.rmssd, color='C3', alpha = 1)
    ax.fill_between(ecg_hr_mean.index, ecg_hr_mean.rmssd-ecg_hr_std.rmssd, ecg_hr_mean.rmssd + ecg_hr_std.rmssd, alpha = 0.3, facecolor='C4')
    plt.gcf().autofmt_xdate()
    ax.set_ylabel("RMSSd [ms]")
    ax.set_xlabel("Date-Time")
    ax.set_xlim([str_time, end_time])
    ax.grid(alpha=0.4)

    #sdnn
    ax = fig.add_subplot(814)
    ax.plot(ecg_hr_mean.index, ecg_hr_mean.sdnn, color='C3', alpha = 1)
    ax.fill_between(ecg_hr_mean.index, ecg_hr_mean.sdnn-ecg_hr_std.sdnn, ecg_hr_mean.sdnn + ecg_hr_std.sdnn, alpha = 0.3, facecolor='C5')
    plt.gcf().autofmt_xdate()
    ax.set_ylabel("SDNN [ms]")
    ax.set_xlabel("Date-Time")
    ax.set_xlim([str_time, end_time])
    ax.grid(alpha=0.4)

    #pnn50
    ax = fig.add_subplot(815)
    ax.plot(ecg_hr_mean.index, ecg_hr_mean.pnn50, color='C3', alpha = 1)
    ax.fill_between(ecg_hr_mean.index, ecg_hr_mean.pnn50-ecg_hr_std.pnn50, ecg_hr_mean.pnn50 + ecg_hr_std.pnn50, alpha = 0.3, facecolor='C6')
    plt.gcf().autofmt_xdate()
    ax.set_ylabel("PNN50 [%]")
    ax.set_xlabel("Date-Time")
    ax.set_xlim([str_time, end_time])
    ax.grid(alpha=0.4)
    
    #Acc Plot
    ax1 = fig.add_subplot(816,sharex=ax)
    
    ax1.plot(acc_data_avg.index, acc_data_avg.accX, color='C1', alpha = 1)
    ax1.fill_between(acc_data_std.index, acc_data_avg.accX-acc_data_std.accX, acc_data_avg.accX + acc_data_std.accX, alpha = 0.3, facecolor='C1')
    plt.gcf().autofmt_xdate()
    ax1.set_ylabel("Acc X [G]")
    ax1.grid(alpha=0.4)

    ax2 = fig.add_subplot(817, sharex=ax)
    ax2.plot(acc_data_avg.index, acc_data_avg.accY, color='C0', alpha = 1)
    ax2.fill_between(acc_data_std.index, acc_data_avg.accY-acc_data_std.accY, acc_data_avg.accY + acc_data_std.accY, alpha = 0.3, facecolor='C0')
    plt.gcf().autofmt_xdate()
    ax2.set_ylabel("Acc Y [G]")
    ax2.grid(alpha=0.4)
    
    ax3 = fig.add_subplot(818, sharex=ax)
    ax3.plot(acc_data_avg.index, acc_data_avg.accZ, color='C2', alpha = 1)
    ax3.fill_between(acc_data_std.index, acc_data_avg.accZ-acc_data_std.accZ, acc_data_avg.accZ+acc_data_std.accZ, alpha = 0.3, facecolor='C2')
    plt.gcf().autofmt_xdate()
    ax3.set_ylabel("Acc Z [G]")
    ax3.grid(alpha=0.4)
    
    plt.show()

def plot_hr_acc_ecg_sz(hr_data_avg, hr_data_std, acc_data_avg, acc_data_std, ecg_hr_mean, ecg_hr_std, sz_start, subject_name, str_time, end_time):
    
    #plot HR for all days
    plt.rcParams.update({'font.size': 18})
    
    fig = plt.figure(figsize=(15,25))
    
    #HR plot
    ax = fig.add_subplot(811)
    ax.plot(hr_data_avg.index, hr_data_avg.heartRate, color='C3', alpha = 1)
    ax.fill_between(hr_data_avg.index, hr_data_avg.heartRate-hr_data_std.heartRate, hr_data_avg.heartRate + hr_data_std.heartRate, alpha = 0.3, facecolor='C3')
    for l in sz_start:
        ax.axvline(x=l, color='C3')
    plt.gcf().autofmt_xdate()
    ax.set_ylabel("HR watch [bpm]")
    ax.set_xlabel("Date-Time")
    ax.set_xlim([str_time, end_time])
    ax.grid(alpha=0.4)
    ax.set_title('Subject : %s' %(subject_name))


    #HR plot ECG
    ax = fig.add_subplot(812)
    ax.plot(ecg_hr_mean.index, ecg_hr_mean.heartRate, color='C3', alpha = 1)
    ax.fill_between(ecg_hr_mean.index, ecg_hr_mean.heartRate-ecg_hr_std.heartRate, ecg_hr_mean.heartRate + ecg_hr_std.heartRate, alpha = 0.3, facecolor='C3')
    for l in sz_start:
        ax.axvline(x=l, color='C3')
    plt.gcf().autofmt_xdate()
    ax.set_ylabel("HR ECG [bpm]")
    ax.set_xlabel("Date-Time")
    ax.set_xlim([str_time, end_time])
    ax.grid(alpha=0.4)

    #rmssd
    ax = fig.add_subplot(813)
    ax.plot(ecg_hr_mean.index, ecg_hr_mean.rmssd, color='C4', alpha = 1)
    ax.fill_between(ecg_hr_mean.index, ecg_hr_mean.rmssd-ecg_hr_std.rmssd, ecg_hr_mean.rmssd + ecg_hr_std.rmssd, alpha = 0.3, facecolor='C4')
    for l in sz_start:
        ax.axvline(x=l, color='C3')
    plt.gcf().autofmt_xdate()
    ax.set_ylabel("RMSSd [ms]")
    ax.set_xlabel("Date-Time")
    ax.set_xlim([str_time, end_time])
    ax.grid(alpha=0.4)

    #sdnn
    ax = fig.add_subplot(814)
    ax.plot(ecg_hr_mean.index, ecg_hr_mean.sdnn, color='C5', alpha = 1)
    ax.fill_between(ecg_hr_mean.index, ecg_hr_mean.sdnn-ecg_hr_std.sdnn, ecg_hr_mean.sdnn + ecg_hr_std.sdnn, alpha = 0.3, facecolor='C5')
    for l in sz_start:
        ax.axvline(x=l, color='C3')
    plt.gcf().autofmt_xdate()
    ax.set_ylabel("SDNN [ms]")
    ax.set_xlabel("Date-Time")
    ax.set_xlim([str_time, end_time])
    ax.grid(alpha=0.4)

    #pnn50
    ax = fig.add_subplot(815)
    ax.plot(ecg_hr_mean.index, ecg_hr_mean.pnn50, color='C6', alpha = 1)
    ax.fill_between(ecg_hr_mean.index, ecg_hr_mean.pnn50-ecg_hr_std.pnn50, ecg_hr_mean.pnn50 + ecg_hr_std.pnn50, alpha = 0.3, facecolor='C6')
    for l in sz_start:
        ax.axvline(x=l, color='C3')
    plt.gcf().autofmt_xdate()
    ax.set_ylabel("PNN50 [%]")
    ax.set_xlabel("Date-Time")
    ax.set_xlim([str_time, end_time])
    ax.grid(alpha=0.4)
    
    #Acc Plot
    ax1 = fig.add_subplot(816,sharex=ax)
    
    ax1.plot(acc_data_avg.index, acc_data_avg.accX, color='C1', alpha = 1)
    ax1.fill_between(acc_data_std.index, acc_data_avg.accX-acc_data_std.accX, acc_data_avg.accX + acc_data_std.accX, alpha = 0.3, facecolor='C1')
    for l in sz_start:
        ax1.axvline(x=l, color='C3')
    plt.gcf().autofmt_xdate()
    ax1.set_ylabel("Acc X [G]")
    ax1.grid(alpha=0.4)

    ax2 = fig.add_subplot(817, sharex=ax)
    ax2.plot(acc_data_avg.index, acc_data_avg.accY, color='C0', alpha = 1)
    ax2.fill_between(acc_data_std.index, acc_data_avg.accY-acc_data_std.accY, acc_data_avg.accY + acc_data_std.accY, alpha = 0.3, facecolor='C0')
    for l in sz_start:
        ax2.axvline(x=l, color='C3')
    plt.gcf().autofmt_xdate()
    ax2.set_ylabel("Acc Y [G]")
    ax2.grid(alpha=0.4)
    
    ax3 = fig.add_subplot(818, sharex=ax)
    ax3.plot(acc_data_avg.index, acc_data_avg.accZ, color='C2', alpha = 1)
    ax3.fill_between(acc_data_std.index, acc_data_avg.accZ-acc_data_std.accZ, acc_data_avg.accZ+acc_data_std.accZ, alpha = 0.3, facecolor='C2')
    for l in sz_start:
        ax3.axvline(x=l, color='C3')
    plt.gcf().autofmt_xdate()
    ax3.set_ylabel("Acc Z [G]")
    ax3.grid(alpha=0.4)
    
    plt.show()



def plot_hr_acc_ecg_circ(hr_data_mean, hr_data_std, acc_data_mean, acc_data_std, ecg_hr_mean, ecg_hr_std, subject_name, str_time, end_time):
    
    
    plt.rcParams.update({'font.size': 12})
    
    fig = plt.figure(figsize=(8,10))
    
    #HR plot
    ax00 = fig.add_subplot(421)
    ax00.plot(hr_data_mean.index, hr_data_mean.heartRate,'o-', color='C3', alpha = 1)
    ax00.fill_between(hr_data_mean.index, hr_data_mean.heartRate-hr_data_std.heartRate, hr_data_mean.heartRate + hr_data_std.heartRate, alpha = 0.3, facecolor='C3')
    ax00.set_ylabel("HR watch\n[bpm]")
    ax00.set_xlabel("Date-Time")
    ax00.grid(alpha=0.4)
    # ax00.set_ylim([40,180 ])


    #HR plot ECG
    ax = fig.add_subplot(422, sharex=ax00)
    ax.plot(ecg_hr_mean.index, ecg_hr_mean.heartRate, 'o-', color='C3', alpha = 1)
    ax.fill_between(ecg_hr_mean.index, ecg_hr_mean.heartRate-ecg_hr_std.heartRate, ecg_hr_mean.heartRate + ecg_hr_std.heartRate, alpha = 0.3, facecolor='C3')
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_ylabel("HR ECG\n[bpm]")
    ax.grid(alpha=0.4)

    #rmssd
    ax = fig.add_subplot(424, sharex=ax00)
    ax.plot(ecg_hr_mean.index, ecg_hr_mean.rmssd,'o-', color='C4', alpha = 1)
    ax.fill_between(ecg_hr_mean.index, ecg_hr_mean.rmssd-ecg_hr_std.rmssd, ecg_hr_mean.rmssd + ecg_hr_std.rmssd, alpha = 0.3, facecolor='C4')
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_ylabel("RMSSd\n[ms]")
    ax.grid(alpha=0.4)

    #sdnn
    ax = fig.add_subplot(426, sharex=ax00)
    ax.plot(ecg_hr_mean.index, ecg_hr_mean.sdnn,'o-', color='C5', alpha = 1)
    ax.fill_between(ecg_hr_mean.index, ecg_hr_mean.sdnn-ecg_hr_std.sdnn, ecg_hr_mean.sdnn + ecg_hr_std.sdnn, alpha = 0.3, facecolor='C5')
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_ylabel("SDNN\n[ms]")
    ax.set_xlabel("Date-Time")
    ax.grid(alpha=0.4)

    #pnn50
    ax = fig.add_subplot(428, sharex=ax00)
    ax.plot(ecg_hr_mean.index, ecg_hr_mean.pnn50,'o-', color='C6', alpha = 1)
    ax.fill_between(ecg_hr_mean.index, ecg_hr_mean.pnn50-ecg_hr_std.pnn50, ecg_hr_mean.pnn50 + ecg_hr_std.pnn50, alpha = 0.3, facecolor='C6')
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_ylabel("PNN50\n[%]")
    ax.set_xlabel("Time of day [hr]")
    ax.grid(alpha=0.4)
    
    #Acc Plot
    ax1 = fig.add_subplot(423, sharex=ax00)
    
    ax1.plot(acc_data_mean.index, acc_data_mean.accX, 'o-', color='C1', alpha = 1)
    ax1.fill_between(acc_data_std.index, acc_data_mean.accX-acc_data_std.accX, acc_data_mean.accX + acc_data_std.accX, alpha = 0.3, facecolor='C1')
    
    # plt.gcf().autofmt_xdate()
    ax1.set_ylabel("Acc X\n[G]")
    ax1.grid(alpha=0.4)

    ax2 = fig.add_subplot(425, sharex=ax00)
    ax2.plot(acc_data_mean.index, acc_data_mean.accY,'o-', color='C0', alpha = 1)
    ax2.fill_between(acc_data_std.index, acc_data_mean.accY-acc_data_std.accY, acc_data_mean.accY + acc_data_std.accY, alpha = 0.3, facecolor='C0')
    
    # plt.gcf().autofmt_xdate()
    ax2.set_ylabel("Acc Y\n[G]")
    ax2.grid(alpha=0.4)
    
    ax3 = fig.add_subplot(427, sharex=ax00)
    ax3.plot(acc_data_mean.index, acc_data_mean.accZ,'o-', color='C2', alpha = 1)
    ax3.fill_between(acc_data_std.index, acc_data_mean.accZ-acc_data_std.accZ, acc_data_mean.accZ+acc_data_std.accZ, alpha = 0.3, facecolor='C2')
    
    # 

    ax3.set_ylabel("Acc Z\n[G]")
    ax3.set_xlabel("Time of day [hr]")
    ax3.grid(alpha=0.4)
    
    sns.despine(left=True, bottom=True, right=True)
    plt.tick_params(top='off', bottom='off', left='off', right='off', labelleft='on', labelbottom='on')
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.1, hspace=0.1)
    plt.gcf().autofmt_xdate()
    plt.suptitle('Subject : %s' %(subject_name))
    plt.show()

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


    #HR plot ECG
    axi = plt.subplot(6,1,2)
    axi.plot(ecg_hr_mean.index, ecg_hr_mean.heartRate, color='C3', alpha = 1)
    axi.fill_between(ecg_hr_mean.index, ecg_hr_mean.heartRate-ecg_hr_std.heartRate, ecg_hr_mean.heartRate + ecg_hr_std.heartRate, alpha = 0.3, facecolor='C3')            
    axi.set_ylabel("HR ECG [bpm]")
    axi.grid(alpha=0.4)

    #rmssd
    axi = plt.subplot(6,1,3)
    axi.plot(ecg_hr_mean.index, ecg_hr_mean.rmssd, color='C4', alpha = 1)
    axi.fill_between(ecg_hr_mean.index, ecg_hr_mean.rmssd-ecg_hr_std.rmssd, ecg_hr_mean.rmssd + ecg_hr_std.rmssd, alpha = 0.3, facecolor='C4')            
    axi.set_ylabel("RMSSD [ms]")
    axi.grid(alpha=0.4)

    #sdnn
    axi = plt.subplot(6,1,4)
    axi.plot(ecg_hr_mean.index, ecg_hr_mean.sdnn, color='C4', alpha = 1)
    axi.fill_between(ecg_hr_mean.index, ecg_hr_mean.sdnn-ecg_hr_std.sdnn, ecg_hr_mean.sdnn + ecg_hr_std.sdnn, alpha = 0.3, facecolor='C4')            
    axi.set_ylabel("SDNN [ms]")
    axi.grid(alpha=0.4)


    #pnn50
    axi = plt.subplot(6,1,5)
    axi.plot(ecg_hr_mean.index, ecg_hr_mean.pnn50, color='C4', alpha = 1)
    axi.fill_between(ecg_hr_mean.index, ecg_hr_mean.pnn50-ecg_hr_std.pnn50, ecg_hr_mean.pnn50 + ecg_hr_std.pnn50, alpha = 0.3, facecolor='C4')            
    axi.set_ylabel("pnn50 [%]")
    axi.grid(alpha=0.4)
    
    #Acc Plot

    #x axis
    axi = plt.subplot(6,1,6)
    axi.plot(acc_data_mean.index, acc_data_mean, color='C0', alpha = 1)
    axi.fill_between(acc_data_mean.index, acc_data_mean-acc_data_std, acc_data_mean + acc_data_std, alpha = 0.3, facecolor='C1')            
    axi.set_ylabel("Acc Mag[G]")
    axi.grid(alpha=0.4)
   

    for axi in ax:
        if sz_start.any():
            for l in sz_start:
                axi.axvline(x=l, color='C3')
        
        ylim = axi.get_ylim()
        if np.any(sleep_data):
            for sl in sleep_data.index:
                axi.fill_betweenx(ylim, sleep_data.loc[sl,'startDate'], sleep_data.loc[sl,'endDate'], alpha=0.2, color ='black')



    plt.xlim([str_time, end_time])
    plt.gcf().autofmt_xdate() 
    
    plt.show()