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


def compute_phase_sync(sig):
    """
    Calculation of phase synchronization using Hilbert transformation
    sensitive to phases, irrespective of the amplitude
    and phase shift, pre-filtering of the signals is necessary

    Parameters
    ----------
    sig: np.array
        2D numpy array of shape (signals, samples), time series ( float)

    Returns
    -------
    phase_sync: float
        ranges between 0 and 1 (1 for the perfect synchronization)

    Example
    -------
    phs = compute_phase_sync(sig)

    References
    ----------
    Quiroga et al. 2008
    """

    if type(sig) != np.ndarray:
        raise TypeError(f"Signals have to be in numpy arrays!")

    # OPTIMIZE: check if this can be done in one array
    sig1_ph = np.angle(hilbert(sig[0]))
    sig2_ph = np.angle(hilbert(sig[1]))

    ph_12 = sig1_ph - sig2_ph
    phase_sync = np.sqrt(np.mean(np.cos(ph_12))**2 + np.mean(np.sin(ph_12))**2)
    # {Quiroga et al. 2008, equation 17 and 18}

    return phase_sync


class PhaseSynchrony(Method):

    algorithm = 'PHASE_SYNCHRONY'
    algorithm_type = 'bivariate'
    version = '1.0.0'
    dtype = [('phase_sync', 'float32')]

    def __init__(self, **kwargs):
        """
        Calculation of phase synchronization using Hilbert transformation
        sensitive to phases, irrespective of the amplitude
        and phase shift, pre-filtering of the signals is necessary

        References
        ----------
        Quiroga et al. 2008
        """
        super().__init__(compute_phase_sync, **kwargs)
