# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.

# Std imports

# Third pary imports
import numpy as np

# Local imports
from ..utils.method import Method
from ..utils.tools import try_jit_decorate


@try_jit_decorate({'nopython': True, 'cache': True})
def compute_signal_stats(sig):
    """
    Function to analyze basic signal statistics

    Parameters:
    ----------
    sig: np.array
        signal to analyze, time series (array, int, float)

    Returns
    -------
    results: tuple
        - power_std: standard deviation of power in band
        - power_mean: mean of power in band
        - power_median: median of power in band
        - power_max: max value of power in band
        - power_min: min value of power in band
        - power_perc25: 25 percentile of power in band
        - power_perc75: 75 percentile of power in band

    Example
    -------
    sig_stats = compute_signal_stats(sig)
    """

    # signal power
    sig_power = sig ** 2

    # compute signal power statistics
    sig_f_pw_std = sig_power.std()
    sig_f_pw_mean = sig_power.mean()
    sig_f_pw_median = np.median(sig_power)
    sig_f_pw_max = sig_power.max()
    sig_f_pw_min = sig_power.min()
    sig_f_pw_perc25 = np.percentile(sig_power, 25)
    sig_f_pw_perc75 = np.percentile(sig_power, 75)

    sig_stats = (sig_f_pw_std, sig_f_pw_mean, sig_f_pw_median, sig_f_pw_max,
                 sig_f_pw_min, sig_f_pw_perc25, sig_f_pw_perc75)

    return sig_stats


class SignalStats(Method):

    algorithm = 'SIGNAL_STATISTICS'
    algorithm_type = 'univariate'
    version = '1.0.0'
    dtype = [('power_std', 'float32'),
             ('power_mean', 'float32'),
             ('power_median', 'float32'),
             ('power_max', 'float32'),
             ('power_min', 'float32'),
             ('power_perc25', 'float32'),
             ('power_perc75', 'float32')]

    def __init__(self, **kwargs):
        """
        Class to analyze basic signal statistics

        Parameters:
        ----------
        sig: np.array
            signal to analyze, time series (array, int, float)

        """

        super().__init__(compute_signal_stats, **kwargs)
