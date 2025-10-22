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


def compute_phase_const(sig, lag=500, lag_step=50):
    """

    **under development**

    calculation of phase consistency between two signals
    irrespective of the amplitude
    pre-filtering of the signals is necessary
    use appropriate lag and step (it calculates phase_const between single
    lag steps in whole length of given time signals)

    Parameters
    ----------
    sig: np.array
        2D numpy array of shape (signals, samples), time series (float)
    lag: int
        negative and positive shift of time series in samples
    lag_step: int
        step of shift in samples

    Returns
    -------
    phase_const: float
        ranges between 0 and 1
        (1 for the phase lock which does not shift during the time period)

    Example
    -------
    phsc = compute_phase_const(sig, 500, 50)
    """

    if type(sig) != np.ndarray:
        raise TypeError(f"Signals have to be in numpy arrays!")

    nstep = int((sig.shape[1] - lag) / lag_step)

    phs_sync_temp = []
    for i in range(0, nstep):
        ind1 = i * lag_step
        ind2 = ind1 + lag

        if ind2 >= sig.shape[1]:
            continue

        # OPTIMIZE: check if we can do this in one array
        sig1_w = sig[0][ind1:ind2]
        sig2_w = sig[1][ind1:ind2]

        sig1_ph = np.unwrap(np.angle(hilbert(sig1_w)))
        sig2_ph = np.unwrap(np.angle(hilbert(sig2_w)))
        ph_12 = sig1_ph - sig2_ph
        phs_sync_temp.append(
            np.sqrt(np.mean(np.cos(ph_12))**2 + np.mean(np.sin(ph_12))**2))

    phase_const = (1 - np.std(phs_sync_temp) / 0.5) * np.mean(phs_sync_temp)

    return phase_const


class PhaseConsistency(Method):

    algorithm = 'PHASE_CONSISTENCY'
    algorithm_type = 'bivariate'
    version = '1.0.0'
    dtype = [('phase_const', 'float32')]

    def __init__(self, **kwargs):
        """
        **under development**

        calculation of phase consistency between two signals
        irrespective of the amplitude
        pre-filtering of the signals is necessary
        use appropriate lag and step (it calculates phase_const between single
        lag steps in whole length of given time signals)

        Parameters
        ----------
        lag: int
            negative and positive shift of time series in samples
        lag_step: int
            step of shift in samples
        """

        super().__init__(compute_phase_const, **kwargs)
