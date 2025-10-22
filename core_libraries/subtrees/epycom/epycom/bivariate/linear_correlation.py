# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.


# Std imports

# Third pary imports
import numpy as np

# Local imports
from ..utils.method import Method


def compute_lincorr(sig, lag=0, lag_step=0):
    """
    Linear correlation (Pearson's coefficient) between two time series

    When lag and lag_step is not 0, shifts the sig[1] from negative
    to positive lag and takes the max correlation (best fit)

    Parameters
    ----------
    sig: np.array
        2D numpy array of shape (signals, samples), time series (int, float)
    lag: int
        negative and positive shift of time series in samples
    lag_step: int
        step of shift

    Returns
    -------
    lincorr: list
        maximum linear correlation in shift
    tau: float
        shift of maximum correlation in samples,
        value in range <-lag,+lag> (float)
        tau<0: sig[1] -> sig[0]
        tau>0: sig[0] -> sig[1]

    Example
    -------
    lincorr,tau = compute_lincorr(sig, 200, 20)
    """

    if type(sig) != np.ndarray:
        raise TypeError(f"Signals have to be in numpy arrays!")

    if lag == 0:
        lag_step = 1
    nstep_lag = int(lag * 2 / lag_step)

    sig1_w = sig[0]
    sig2_w = sig[1]

    sig1_wl = sig1_w[lag:len(sig1_w) - lag]

    lincorr = []
    for i in range(0, nstep_lag + 1):
        ind1 = i * lag_step
        ind2 = ind1 + len(sig1_wl)

        sig2_wl = sig2_w[ind1:ind2]

        corr_val = np.corrcoef(sig1_wl, sig2_wl)
        lincorr.append(corr_val[0][1])

    return np.max(lincorr), lincorr.index(max(lincorr))


class LinearCorrelation(Method):

    algorithm = 'LINEAR_CORRELATION'
    algorithm_type = 'bivariate'
    version = '1.0.0'
    dtype = [('max_corr', 'float32'),
             ('tau', 'float32')]

    def __init__(self, **kwargs):
        """
        Linear correlation (Pearson's coefficient) between two time series

        When win and win_step is not 0, calculates evolution of correlation

        When win>len(sig) or win<=0, calculates only one corr coef

        When lag and lag_step is not 0, shifts the sig[1] from negative
        to positive lag and takes the max correlation (best fit)

        Parameters
        ----------
        lag: int
            negative and positive shift of time series in samples
        lag_step: int
            step of shift

        """
        super().__init__(compute_lincorr, **kwargs)
