# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering;
# Institute of Scientific Instruments of the CAS, v. v. i., Medical signals -
# Computational neuroscience. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.

# Std imports

# Third pary imports
import numpy as np
from scipy import signal

# Local imports
from ..utils.method import Method
#from ..utils.tools import try_jit_decorate


def compute_powerline_noise(sig, fs, freq=60):

    
    """
    Function to detect the proportio of line noise in the signal

    Parameters:
    ----------
    sig: np.array
        signal to analyze, time series (array, int, float)
    fs: float
        sampling frequency
    freq: float or int
        line frequency (50 or 60 Hz)

    Returns
    -------
    score: float
        line frequency score

    Example
    -------
    sat = compute_saturation(sig)
    """

    #low pass pre filter
    b, a = signal.butter(3, 5/(fs/2), btype='high',analog=False)
    x = signal.filtfilt(b, a, sig)
    
    # signal normalization
    x = (x - np.mean(x))/np.std(x)
    
    
    # notch filter 50/60 Hz, bandwidth 2Hz
    b, a = signal.iirnotch(freq, freq/1, fs)
    x_filt = signal.filtfilt(b, a, x)
    
    subtraction = np.abs(x - x_filt)
    
    
    return np.std(subtraction)/np.std(x)


class PowerlineNoise(Method):

    algorithm = 'POWERLINE_NOISE'
    algorithm_type = 'artifact'
    version = '1.0.0'
    dtype = [('score', 'float32')]

    def __init__(self, **kwargs):
        """
        Class to detect flat line (saturation or missing data)

        Parameters:
        ----------
        sig: np.array
            signal to analyze, time series (array, int, float)

        """

        super().__init__(compute_powerline_noise, **kwargs)
