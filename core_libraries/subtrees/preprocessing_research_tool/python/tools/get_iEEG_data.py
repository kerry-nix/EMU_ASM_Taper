# pylint: disable-msg=C0103
import ieeg
from ieeg.auth import Session
import pandas as pd
import pickle


# from .pull_patient_localization import pull_patient_localization
# from pull_patient_localization import pull_patient_localization
from numbers import Number
import numpy as np
import time, os, warnings
from beartype import beartype
from typing import Union

from .clean_labels import clean_labels


def _pull_iEEG(ds: ieeg.dataset.Dataset, start_usec: Number, duration_usec: Number, channel_ids: list) -> np.ndarray:
    """
    Pull data while handling iEEGConnectionError
    """
    while True:
        try:
            data = ds.get_data(start_usec, duration_usec, channel_ids)
            return data
        except Exception as e:
            if '500' in str(e) | '502' in str(e) | '503' in str(e) | '504' in str(e):
                time.sleep(1)
            else:
                raise e


@beartype
def get_ieeg_data(
    username: str,
    password_bin_file: str,
    iEEG_filename: str,
    start_time: Number,
    stop_time: Number,
    select_elecs: list[Union[str, int]] = None,
    ignore_elecs: list[Union[str, int]] = None,
    outputfile=None,
):
    """ "
    2020.04.06. Python 3.7
    Andy Revell, adapted by Akash Pattnaik (2021.06.23)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Purpose:
    To get iEEG data from iEEG.org. Note, you must download iEEG python package from GitHub - instructions are below
    1. Gets time series data and sampling frequency information. Specified electrodes are removed.
    2. Saves as a pickle format
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Input
        username: your iEEG.org username
        password_bin_file: your iEEG.org password bin_file
        iEEG_filename: The file name on iEEG.org you want to download from
        start_time: the start time in the iEEG_filename. In seconds
        stop_time: the stop time in the iEEG_filename. In seconds.
            iEEG.org needs a duration input: this is calculated by stop_time - start_time
        select_elecs: the electrode/channel names/indices you want to select.
        ignore_elecs: the electrode/channel names/indices you want to exclude. EXACT MATCH on iEEG.org. Caution: some may be LA08 or LA8
        outputfile: the path and filename you want to save.
            PLEASE INCLUDE EXTENSION .pickle.
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Output:
        Saves file outputfile as a pickle. For more info on pickling, see https://docs.python.org/3/library/pickle.html
        Briefly: it is a way to save + compress data. it is useful for saving lists, as in a list of time series data and sampling frequency together along with channel names
        List index 0: Pandas dataframe. T x C (rows x columns). T is time. C is channels.
        List index 1: float. Sampling frequency. Single number
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Example usage:
    username = 'arevell'
    password = 'password'
    iEEG_filename='HUP138_phaseII'
    start_time = 248432.34
    stop_time = 248525.74
    removed_channels = ['EKG1', 'EKG2', 'CZ', 'C3', 'C4', 'F3', 'F7', 'FZ', 'F4', 'F8', 'LF04', 'RC03', 'RE07', 'RC05', 'RF01', 'RF03', 'RB07', 'RG03', 'RF11', 'RF12']
    outputfile = '/Users/andyrevell/mount/DATA/Human_Data/BIDS_processed/sub-RID0278/eeg/sub-RID0278_HUP138_phaseII_248432340000_248525740000_EEG.pickle'
    get_iEEG_data(username, password, iEEG_filename, start_time_usec, stop_time_usec, removed_channels, outputfile)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    To run from command line:
    python3.6 -c 'import get_iEEG_data; get_iEEG_data.get_iEEG_data("arevell", "password", "HUP138_phaseII", 248432340000, 248525740000, ["EKG1", "EKG2", "CZ", "C3", "C4", "F3", "F7", "FZ", "F4", "F8", "LF04", "RC03", "RE07", "RC05", "RF01", "RF03", "RB07", "RG03", "RF11", "RF12"], "/gdrive/public/DATA/Human_Data/BIDS_processed/sub-RID0278/eeg/sub-RID0278_HUP138_phaseII_D01_248432340000_248525740000_EEG.pickle")'
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #How to get back pickled files
    with open(outputfile, 'rb') as f: data, fs = pickle.load(f)
    """

    # print("\n\nGetting data from iEEG.org:")
    # print("iEEG_filename: {0}".format(iEEG_filename))
    # print("start_time_usec: {0}".format(start_time_usec))
    # print("stop_time_usec: {0}".format(stop_time_usec))
    # print("ignore_elecs: {0}".format(ignore_elecs))
    # if outputfile:
    #     print("Saving to: {0}".format(outputfile))
    # else:
    #     print("Not saving, returning data and sampling frequency")

    # Pull and format metadata from patient_localization_mat

    # Added by Haoer
    current_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    password_bin_file = os.path.join(current_dir,password_bin_file)
    pwd = open(password_bin_file, "r").read()

    assert start_time < stop_time, "CNTtools:invalidTimeRange"
    assert start_time >= 0, "CNTtools:invalidTimeRange"
    start_time_usec = int(start_time * 1e6)
    stop_time_usec = int(stop_time * 1e6)
    duration = stop_time_usec - start_time_usec

    # check
    # no meta data checking currently
    #if os.path.exists("meta_data.csv"):
    #    meta_data = pd.read_csv("meta_data.csv")
    #    assert iEEG_filename in meta_data["filename"].values, "CNTtools:invalidFileName"
    #    assert (
    #        duration <= meta_data["duration"][meta_data["filename"].eq(iEEG_filename)]
    #    ), "CNTtools:invalidTimeRange"
    #else:
    while True:
        try:
            s = Session(username, pwd)
            ds = s.open_dataset(iEEG_filename)
            all_channel_labels = ds.get_channel_labels()
            break
        except Exception as e:
            if 'Authentication' in str(e):
                raise AssertionError("CNTtools:invalidLoginInfo")
            elif '404' in str(e) or 'NoSuchDataSnapshot' in str(e):
                raise AssertionError("CNTtools:invalidFileName")
            elif '500' in str(e) | '502' in str(e) | '503' in str(e) | '504' in str(e):
                time.sleep(1)
            else:
                raise e
    assert len(all_channel_labels) > 0, "CNTtools:emptyFile"
    end_sec = ds.get_time_series_details(all_channel_labels[0]).duration
    assert stop_time_usec <= end_sec, "CNTtools:invalidTimeRange"       
    all_channel_labels = clean_labels(all_channel_labels)

    if select_elecs is not None:
        elec_type = type(select_elecs[0])
        assert all(isinstance(i, elec_type) for i in select_elecs), "CNTtools:invalidElectrodeList"
        if elec_type == int:
            channel_ids = [i for i in select_elecs if i >= 0 & i < len(all_channel_labels)]
            if len(channel_ids) < len(select_elecs):
                warnings.warn("CNTtools:invalidChannelID, invalid channels ignored.")
            channel_names = [all_channel_labels[e] for e in channel_ids]
        elif elec_type == str:
            select_elecs = clean_labels(select_elecs)
            channel_ids = [
                i for i, e in enumerate(all_channel_labels) if e in select_elecs
            ]
            if len(channel_ids) < len(select_elecs):
                warnings.warn("CNTtools:invalidChannelID, invalid channels ignored.")
            channel_names = [all_channel_labels[e] for e in channel_ids]
        else:
            print("Electrodes not given as a list of ints or strings")

    elif ignore_elecs is not None:
        elec_type = type(ignore_elecs[0])
        assert all(isinstance(i, elec_type) for i in ignore_elecs), "CNTtools:invalidElectrodeList"
        if elec_type == int:
            channel_ids = [
                i
                for i in np.arange(len(all_channel_labels))
                if i not in ignore_elecs
            ]
            if len(channel_ids) > len(all_channel_labels) - len(ignore_elecs):
                warnings.warn("CNTtools:invalidChannelID, invalid channels ignored.")
            channel_names = [all_channel_labels[e] for e in channel_ids]
        elif elec_type == str:
            ignore_elecs = clean_labels(ignore_elecs)
            channel_ids = [
                i
                for i, e in enumerate(all_channel_labels)
                if e not in ignore_elecs
            ]
            if len(channel_ids) > len(all_channel_labels) - len(ignore_elecs):
                warnings.warn("CNTtools:invalidChannelID, invalid channels ignored.")
            channel_names = [
                e for e in all_channel_labels if e not in ignore_elecs
            ]
        else:
            print("Electrodes not given as a list of ints or strings")

    else:
        channel_ids = np.arange(len(all_channel_labels))
        channel_names = all_channel_labels

    try:
        data = ds.get_data(start_time_usec, duration, channel_ids)
    except Exception as e:
        # clip is probably too big, pull chunks and concatenate
        clip_size = 60 * 1e6
        clip_start = start_time_usec
        data = None
        while clip_start + clip_size < stop_time_usec:
            if data is None:
                # data = ds.get_data(clip_start, clip_size, channel_ids)
                data = _pull_iEEG(ds, clip_start, clip_size, channel_ids)
            else:
                # data = np.concatenate(([data, ds.get_data(clip_start, clip_size, channel_ids)]), axis=0)
                data = np.concatenate(
                    ([data, _pull_iEEG(ds, clip_start, clip_size, channel_ids)]), axis=0
                )
            clip_start = clip_start + clip_size
        # data = np.concatenate(([data, ds.get_data(clip_start, stop_time_usec - clip_start, channel_ids)]), axis=0)

    # df = pd.DataFrame(data, columns=channel_names)
    fs = ds.get_time_series_details(ds.ch_labels[0]).sample_rate  # get sample rate

    if outputfile:
        with open(outputfile, "wb") as f:
            pickle.dump([data, channel_names, fs], f)
    else:
        return data, np.array(channel_names), fs

    # session.delete
    # clear variables


""""
Download and install iEEG python package - ieegpy
GitHub repository: https://github.com/ieeg-portal/ieegpy
If you downloaded this code from https://github.com/andyrevell/paper001.git then skip to step 2
1. Download/clone ieepy. 
    git clone https://github.com/ieeg-portal/ieegpy.git
2. Change directory to the GitHub repo
3. Install libraries to your python Path. If you are using a virtual environment (ex. conda), make sure you are in it
    a. Run:
        python setup.py build
    b. Run: 
        python setup.py install
              
"""
