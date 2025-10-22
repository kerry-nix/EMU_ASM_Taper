# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.

# Std imports

# Third pary imports
import numpy as np
import pandas as pd
from scipy import stats as stats

# Local imports
from ..utils.method import Method


def compute_shanon_entropy(sig):
    """
    Fucntion computes shannon entropy of given signal

    Parameters
    ----------
    sig: np.ndarray
        Signal to analyze

    Returns
    -------
    entro: np.float64
        Computed Shannon entropy of given signal
    """
    pd_series = pd.Series(sig)
    counts = pd_series.value_counts()
    entro = stats.entropy(counts, base=2)               # shan_en = -sum(p(xi)*log(p(xi)))
    return entro


class ShannonEntropy(Method):

    algorithm = 'SHANNON_ENTROPY'
    algorithm_type = 'univariate'
    version = '1.0.0'
    dtype = [('shannon', 'float32')]

    def __init__(self, **kwargs):
        """
        Shannon entropy

        Parameters
        ----------
        sig: np.ndarray
            Signal to analyze
        """

        super().__init__(compute_shanon_entropy, **kwargs)
        self._event_flag = False

