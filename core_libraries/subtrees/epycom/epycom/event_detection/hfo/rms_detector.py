# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.


# Std imports

# Third pary imports
import numpy as np

# Local imports
from ...utils.signal_transforms import compute_rms
from ...utils.thresholds import th_std
from ...utils.method import Method


def detect_hfo_rms(sig, fs=5000, threshold=3, window_size=100,
                   window_overlap=0.25):
    """
    Root mean square detection algorithm {Staba et al. 2002,
    Blanco et al 2010}.

    Parameters
    ----------
    sig: np.ndarray
        1D array with raw data (already filtered if required)
    fs: int
        Sampling frequency
    threshold: float
        Number of standard deviations to use as a threshold
    window_size: int
        Sliding window size in samples
    window_overlap: float
        Fraction of the window overlap (0 to 1)

    Returns
    -------
    output: list
        List of tuples with the following structure:
        (event_start, event_stop)

    References
    ----------
    [1] R. J. Staba, C. L. Wilson, A. Bragin, I. Fried, and J. Engel,
    “Quantitative Analysis of High-Frequency Oscillations (80 − 500 Hz)
    Recorded in Human Epileptic Hippocampus and Entorhinal Cortex,”
    J. Neurophysiol., vol. 88, pp. 1743–1752, 2002.
    """

    # Calculate window values for easier operation
    window_increment = int(np.ceil(window_size * window_overlap))

    output = []

    # Overlapping window

    win_start = 0
    win_stop = window_size
    n_windows = int(np.ceil((len(sig) - window_size) / window_increment)) + 1
    RMS = np.empty(n_windows)
    RMS_i = 0
    while win_start < len(sig):
        if win_stop > len(sig):
            win_stop = len(sig)

        RMS[RMS_i] = compute_rms(sig[int(win_start):int(win_stop)],
                                 window_size)[0]

        if win_stop == len(sig):
            break

        win_start += window_increment
        win_stop += window_increment

        RMS_i += 1

    # Create threshold
    det_th = th_std(RMS, threshold)

    # Detect
    RMS_idx = 0
    while RMS_idx < len(RMS):
        if RMS[RMS_idx] >= det_th:
            event_start = RMS_idx * window_increment
            while RMS_idx < len(RMS) and RMS[RMS_idx] >= det_th:
                RMS_idx += 1
            event_stop = (RMS_idx * window_increment) + window_size

            if event_stop > len(sig):
                event_stop = len(sig)

            # Optional feature calculations can go here

            # Write into output
            output.append((event_start, event_stop))

            RMS_idx += 1
        else:
            RMS_idx += 1

    return output


class RootMeanSquareDetector(Method):

    algorithm = 'ROOTMEANSQUARE_DETECTOR'
    algorithm_type = 'event'
    version = '1.0.0'
    dtype = [('event_start', 'int32'),
             ('event_stop', 'int32')]

    def __init__(self, **kwargs):
        """
        Root mean square detection algorithm.

        Parameters
        ----------
        fs: int
            Sampling frequency
        threshold: float
            Number of standard deviations to use as a threshold
        window_size: int
            Sliding window size in samples
        window_overlap: float
            Fraction of the window overlap (0 to 1)
        sample_offset: int
            Offset which is added to the final detection. This is used when the
            function is run in separate windows. Default = 0

        References
        ----------
        [1] R. J. Staba, C. L. Wilson, A. Bragin, I. Fried, and J. Engel,
        “Quantitative Analysis of High-Frequency Oscillations (80 − 500 Hz)
        Recorded in Human Epileptic Hippocampus and Entorhinal Cortex,”
        J. Neurophysiol., vol. 88, pp. 1743–1752, 2002.
        """

        super().__init__(detect_hfo_rms, **kwargs)
