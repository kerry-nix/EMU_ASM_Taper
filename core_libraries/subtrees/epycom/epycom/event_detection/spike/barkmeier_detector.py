# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.


# Std imports

# Third pary imports
import numpy as np
from scipy.signal import butter, filtfilt, find_peaks

# Local imports
from ...utils.method import Method


def detect_spikes_barkmeier(sig, fs=5000, scale=70, std_coeff=4,
                            through_search=0.05,
                            det_thresholds={'LS': 700,
                                            'RS': 700,
                                            'TAMP': 600,
                                            'LD': 0.01,
                                            'RD': 0.01},
                            filter_spec={'narrow': [20, 50],
                                         'broad': [1, 80]},
                            win_idx=None):
    """
    Python version of Barkmeier's EEG spike detector. {Barkmeier et al. 2011}

    Parameters
    ----------
    sig: np.ndarray
        1D numpy array of EEG data
    fs: int
        sampling frequency of the signal
    scale: float\
        scaling parameter (def=70)
    std_coeff: float
        z-score threshold for spike detection (def=4)
    through_search: float
        extent to which search for spike throughs in s (def=0.04)
    det_thresholds: dict
        detection thresholds (dictionary)
        {'LS':700, # Left slope
         'RS':700, # Right slope
         'TAMP':600, # Total amplitude
         'LD':0.01, # Left duration
         'RD':0.01} # Right duration
    filter_spec: dict
        narrow and broad band filter specifications
        {'narrow':[20, 50],
         'broad':[1, 80]}
    win_idx: int
        Statistical window index. This is used when the
        function is run in separate windows. Default = None

    Returns
    -------
    output: list
        List of tuples with the following structure of detections:
        (event_peak, event_amp, left_amp, left_dur, right_amp, right_dur)
    """

    # Create filter coeficients

    bh1, ah1 = butter(2, filter_spec['narrow'][0] / (fs / 2), 'highpass')
    bl1, al1 = butter(4, filter_spec['narrow'][1] / (fs / 2), 'lowpass')
    bh2, ah2 = butter(2, filter_spec['broad'][0] / (fs / 2),  'highpass')
    bl2, al2 = butter(4, filter_spec['broad'][1] / (fs / 2), 'lowpass')

    output = []

    last_idx = -0.005 * fs

    # Filter data
    fx_narrow = filtfilt(bh1, ah1, sig)
    fx_narrow = filtfilt(bl1, al1, fx_narrow)
    fx_broad = filtfilt(bh2, ah2, sig)
    fx_broad = filtfilt(bl2, al2, fx_broad)

    # Scale the data
    scale_factor = scale / np.median(np.mean(np.abs(fx_broad)))
    fx_broad *= scale_factor

    thresh = np.mean(np.abs(fx_narrow)) + std_coeff * np.std(np.abs(fx_narrow))
    peak_idxs = np.where(fx_narrow > thresh)[0]
    peaks = fx_narrow[peak_idxs]

    pis = peak_idxs[find_peaks(peaks)[0]]  # Getting the maxima

    # Run through peaks and calculate slopes and threshold them
    for pi in pis:

        # Get correct spike index and voltage
        l_idx = int(pi - fs * 0.002)
        r_idx = int(pi + fs * 0.002)
        if l_idx < 0:
            l_idx = 0
        if r_idx > len(sig):
            r_idx = len(sig)

        spike_i = np.argmax(fx_broad[l_idx:r_idx])
        spike_i += l_idx
        spike_V = fx_broad[spike_i]

        # Get the left trough index and voltage
        l_idx = spike_i - int(fs * through_search)
        if l_idx < 0:
            l_idx = 0
        if spike_i == l_idx:
            continue
        left_i = np.argmin(fx_broad[l_idx:spike_i])
        left_i += l_idx
        left_V = fx_broad[left_i]

        # Get the right through index and voltage
        r_idx = spike_i + int(fs * through_search)
        if r_idx < 0:
            r_idx = len(sig)
        if spike_i == r_idx:
            continue
        right_i = np.argmin(fx_broad[spike_i:r_idx])
        right_i += spike_i
        right_V = fx_broad[right_i]

        # Get amp, dur and slope of the left halfwave
        l_amp = spike_V - left_V
        l_dur = (spike_i - left_i) / fs
        l_slope = l_amp / l_dur

        # Get amp, dur and slope of the right halfwave
        r_amp = spike_V - right_V
        r_dur = (right_i - spike_i) / fs
        r_slope = r_amp / r_dur

        # Threshold
        if (((l_slope > det_thresholds['LS'] and
              r_slope > det_thresholds['RS'] and
              l_amp + r_amp > det_thresholds['TAMP'] and
              l_dur > det_thresholds['LD'] and
              r_dur > det_thresholds['RD'])
             or
             (l_slope < det_thresholds['LS'] and
              r_slope < det_thresholds['RS'] and
              l_amp + r_amp < det_thresholds['TAMP'] and
              l_dur > det_thresholds['LD'] and
              r_dur > det_thresholds['RD']))
                and spike_i - last_idx > 0.005):
            if win_idx is not None:
                output.append((int(spike_i), spike_V,
                               l_amp, l_dur,
                               r_amp, r_dur,
                               win_idx))
            else:
                output.append((int(spike_i), spike_V,
                               l_amp, l_dur,
                               r_amp, r_dur))
            last_idx = spike_i

    return output


class BarkmeierDetector(Method):

    algorithm = 'BARKMEIER_DETECTOR'
    algorithm_type = 'event'
    version = '1.0.0'
    dtype = [('event_peak', 'int32'),
             ('event_amp', 'float32'),
             ('left_amp', 'float32'),
             ('left_dur', 'float32'),
             ('right_amp', 'float32'),
             ('right_dur', 'float32')]

    def __init__(self, **kwargs):
        """
        Python version of Barkmeier's EEG spike detector.
         {Barkmeier et al. 2011}

        Parameters
        ----------
        fs: int
            sampling frequency of the signal
        scale: float\
            scaling parameter (def=70)
        std_coef: float
            z-score threshold for spike detection (def=4)
        through_search: float
            extent to which search for spike throughs in s (def=0.04)
        det_thresholds: dict
            detection thresholds (dictionary)
            {'LS':700, # Left slope
             'RS':700, # Right slope
             'TAMP':600, # Total amplitude
             'LD':0.01, # Left duration
             'RD':0.01} # Right duration
        filter_spec: dict
            narrow and broad band filter specifications
            {'narrow':[20, 50],
             'broad':[1, 80]}
        sample_offset: int
            Offset which is added to the final detection. This is used when the
            function is run in separate windows. Default = 0
        """

        super().__init__(detect_spikes_barkmeier, **kwargs)
