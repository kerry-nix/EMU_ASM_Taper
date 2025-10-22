# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.

# Std imports
from math import isclose

# Third pary imports
from scipy.signal import butter, filtfilt

# Local imports
from epycom.event_detection import BarkmeierDetector

from epycom.event_detection import (LineLengthDetector,
                                    RootMeanSquareDetector,
                                    HilbertDetector,
                                    CSDetector)


# ----- Spikes -----
def test_detect_spikes(create_testing_eeg_data, benchmark):
    compute_instance = BarkmeierDetector()
    dets = benchmark(compute_instance.run_windowed,
                     create_testing_eeg_data, 50000)
    compute_instance.run_windowed(create_testing_eeg_data,
                                  5000,
                                  n_cores=2)

    expected_vals = (20242,
                     1368.2334,
                     1517.9938,
                     0.05,
                     1486.8751,
                     0.0376)

    for exp_val, det in zip(expected_vals, dets[0]):
        assert isclose(det, exp_val, abs_tol=10e-5)


# ----- HFO -----
def test_detect_hfo_ll(create_testing_eeg_data, benchmark):
    fs = 5000
    b, a = butter(3, [80 / (fs / 2), 600 / (fs / 2)], 'bandpass')
    filt_data = filtfilt(b, a, create_testing_eeg_data)
    window_size = int((1 / 80) * fs)

    compute_instance = LineLengthDetector()
    compute_instance.params = {'window_size': window_size}
    dets = benchmark(compute_instance.run_windowed,
                     filt_data, 50000)

    compute_instance.run_windowed(filt_data,
                                  5000,
                                  n_cores=2)

    expected_vals = [(5040, 5198),
                     (34992, 35134)]

    for exp_val, det in zip(expected_vals, dets):
        assert det[0] == exp_val[0]
        assert det[1] == exp_val[1]


def test_detect_hfo_rms(create_testing_eeg_data, benchmark):
    fs = 5000
    b, a = butter(3, [80 / (fs / 2), 600 / (fs / 2)], 'bandpass')
    filt_data = filtfilt(b, a, create_testing_eeg_data)
    window_size = int((1 / 80) * fs)

    compute_instance = RootMeanSquareDetector()
    compute_instance.params = {'window_size': window_size}
    dets = benchmark(compute_instance.run_windowed,
                     filt_data, 50000)

    compute_instance.run_windowed(filt_data,
                                  5000,
                                  n_cores=2)

    expected_vals = [(5040, 5198),
                     (35008, 35134)]

    for exp_val, det in zip(expected_vals, dets):
        assert det[0] == exp_val[0]
        assert det[1] == exp_val[1]


def test_detect_hfo_hilbert(create_testing_eeg_data, benchmark):
    compute_instance = HilbertDetector()
    compute_instance.params = {'fs': 5000,
                               'low_fc': 80,
                               'high_fc': 600,
                               'threshold': 7}
    dets = benchmark(compute_instance.run_windowed,
                     create_testing_eeg_data, 50000)

    compute_instance.run_windowed(create_testing_eeg_data,
                                  5000,
                                  n_cores=2)

    expected_vals = [(5056, 5123),
                     (35028, 35063)]

    for exp_val, det in zip(expected_vals, dets):
        assert det[0] == exp_val[0]
        assert det[1] == exp_val[1]


def test_detect_hfo_cs_beta(create_testing_eeg_data, benchmark):
    compute_instance = CSDetector()
    compute_instance.params = {'fs': 5000,
                               'low_fc': 40,
                               'high_fc': 1000,
                               'threshold': 0.1,
                               'cycs_per_detect': 4.0}

    dets = benchmark(compute_instance.run_windowed,
                     create_testing_eeg_data, 50000)

    compute_instance.run_windowed(create_testing_eeg_data,
                                  5000,
                                  n_cores=2)

    # Only the second HFO is caught by CS (due to signal artificiality)
    expected_vals = [(34992, 35090),  # Band detection
                     (34992, 35090)]  # Conglomerate detection

    for exp_val, det in zip(expected_vals, dets):
        assert det[0] == exp_val[0]
        assert det[1] == exp_val[1]
