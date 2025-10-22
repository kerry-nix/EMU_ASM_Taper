# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.

# Third pary imports
import numpy as np
from numba import njit

# Local imports
from ..utils.method import Method


@njit('f8(f8[:], f8[:])', cache=True)
def _maxdist(x_i, x_j):
    dist = 0

    leni = len(x_i)
    lenj = len(x_j)

    if leni < lenj:
        n = len(x_i)
    else:
        n = len(x_j)

    for ua in range(n):
        if abs(x_i[ua] - x_j[ua]) > dist:
            dist = abs(x_i[ua] - x_j[ua])

    return dist


@njit('f8(i8, i8, f8, f8[:])', cache=True)
def _phi_jitted(m, N, r, sig):
    z = N - m + 1

    xlen = N - m + 1
    x = np.full((xlen, m), np.inf, dtype='float64')

    # Sampling the signal
    for i in range(xlen):
        x[i] = sig[i: i + m]

    C = np.full(len(sig), np.inf, dtype='float64')
    iterator = cnt = 0
    for x_i in x:
        for x_j in x:
            if _maxdist(x_i, x_j) <= r:
                cnt += 1
        C[iterator] = cnt / (N - m + 1.0)
        cnt = 0
        iterator += 1

    C = C[:iterator]

    phi = 0
    for c in C:
        phi = phi+np.log(c)

    return phi/z


@njit('f8(f8[:], f8, i8)', cache=True)
def compute_approximate_entropy(sig, r, m):
    """
    Function computes approximate entropy of given signal

    Parameters
    ----------
    sig: np.ndarray
        1D signal
    r: np.float64
        filtering treshold, recommended values: (0.1-0.25)*np.nanstd(sig)
    m: int
        window length of compared run of data, recommended (2-8)

    Returns
    -------
    entro: numpy.float64
        approximate entropy

    Example
    -------
    signal_entropy = approximate_entropy(data, 0.1*np.nanstd(data))
    """

    N = sig.shape[0]
    return abs(_phi_jitted(m + 1, N, r, sig) - _phi_jitted(m, N, r, sig))


class ApproximateEntropy(Method):

    algorithm = 'APPROXIMATE_ENTROPY'
    algorithm_type = 'univariate'
    version = '1.0.0'
    dtype = [('apen', 'float32')]

    def __init__(self, **kwargs):
        """
        Approximate entropy

        Parameters
        ----------
        sig: np.ndarray
            1D signal
        m: int
            window length of compared run of data, recommended (2-8)
        r: float64
            filtering treshold, recommended values: (0.1-0.25)*std
       """

        super().__init__(compute_approximate_entropy, **kwargs)
        self._event_flag = False
