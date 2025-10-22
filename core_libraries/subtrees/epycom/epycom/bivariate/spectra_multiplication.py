# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.


# Std imports

# Third pary imports
import numpy as np
from scipy.signal import hilbert

# Local imports
from ..utils.method import Method


def compute_spect_multp(sig):
    """
    Multiply spectra of two time series and transforms it back to time domain,
    where the mean and std is calculated

    Parameters
    ----------
    sig: np.array
        2D numpy array of shape (signals, samples), time series (int, float)

    Returns
    -------
    sig_sm_mean: float
        aritmetic mean value of multiplied signals
    sig_sm_std: float
        standard deviation of multiplied signals

    Example
    -------
    mspect = spect_multp(sig)
    """

    if type(sig) != np.ndarray:
        raise TypeError(f"Signals have to be in numpy arrays!")

    # OPTIMIZE: check if we can do this in 1 array instead of 2
    fft_1 = np.fft.rfft(sig[0])
    fft_2 = np.fft.rfft(sig[1])
    fft_p = np.multiply(fft_1, fft_2)

    sig_sm = np.abs(hilbert(np.fft.irfft(fft_p)))

    sig_sm_mean = np.mean(sig_sm)
    sig_sm_std = np.std(sig_sm)

    return sig_sm_mean, sig_sm_std


class SpectraMultiplication(Method):

    algorithm = 'SPECTRA_MULTIPLICATION'
    algorithm_type = 'bivariate'
    version = '1.0.0'
    dtype = [('sm_mean', 'float32'),
             ('sm_std', 'float32')]

    def __init__(self, **kwargs):
        """
        Multiply spectra of two time series and transforms it back to time
        domain where the mean and std is calculated
        """
        super().__init__(compute_spect_multp, **kwargs)
