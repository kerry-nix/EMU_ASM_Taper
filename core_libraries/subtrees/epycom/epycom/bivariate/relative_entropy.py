# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.


# Std imports

# Third pary imports
import numpy as np
from scipy.stats import entropy

# Local imports
from ..utils.method import Method


def compute_relative_entropy(sig):
    """
    Calculation of Kullback-Leibler divergence:
    relative entropy of sig[0] with respect to sig[1]
    and relative entropy of sig[1] with respect to sig[0]

    Parameters
    ----------
    sig: np.array
        2D numpy array of shape (signals, samples), time series (int, float)

    Returns
    -------
    ren: float
        max value of relative entropy between sig[0] and sig[1]

    Example:
    -------
    ren = compute_relative_entropy(sig)
    """

    if type(sig) != np.ndarray:
        raise TypeError(f"Signals have to be in numpy arrays!")

    # OPTIMIZE - check if we can do this in one array
    h1 = np.histogram(sig[0], 10)
    h2 = np.histogram(sig[1], 10)

    ren = entropy(h1[0], h2[0])
    ren21 = entropy(h2[0], h1[0])

    if ren21 > ren:
        ren = ren21

    if ren == float('Inf'):
        ren = np.nan

    return ren


class RelativeEntropy(Method):

    algorithm = 'RELATIVE_ENTROPY'
    algorithm_type = 'bivariate'
    version = '1.0.0'
    dtype = [('ren', 'float32')]

    def __init__(self, **kwargs):
        """
        Calculation of Kullback-Leibler divergence:
        relative entropy of sig1 with respect to sig2
        and relative entropy of sig2 with respect to sig1
        """

        super().__init__(compute_relative_entropy, **kwargs)
