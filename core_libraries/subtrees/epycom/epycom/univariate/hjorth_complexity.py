# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.

# Std imports

# Third pary imports
import numpy as np

# Local imports
from .hjorth_mobility import compute_hjorth_mobility
from ..utils.method import Method


def compute_hjorth_complexity(signal, fs=5000):
    """
    Compute Hjorth complexity of time series

    Parameters
    ----------
    signal: np.array
        Signal to analyze, time series (array, int, float)
    fs: float
        Sampling frequency of the time series

    Returns
    -------
    hjorth_complexity: float

    Example
    -------
    hjorth_complexity = compute_hjorth_complexity(data, 5000)

    Note
    ----
    result is NOT frequency dependent
    """

    mob = compute_hjorth_mobility(signal, fs)
    # diff signal is one sample shorter
    mobd = compute_hjorth_mobility(np.diff(signal) * fs, fs)
    hjorth_complexity = mobd / mob
    return hjorth_complexity


class HjorthComplexity(Method):

    algorithm = 'HJORTH_COMPLEXITY'
    algorithm_type = 'univariate'
    ersion = '1.0.0'
    dtype = [('hjorth_complexity', 'float32')]

    def __init__(self, **kwargs):
        """
        Hjorth complexity of time series

        Parameters
        ----------
        fs: float
            Sampling frequency of the time series

        Note
        ----
        result is NOT frequency dependent
        """

        super().__init__(compute_hjorth_complexity, **kwargs)
