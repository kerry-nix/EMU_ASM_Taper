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


@njit('f8(f8[:], f8, i8)', cache=True)
def compute_sample_entropy(sig, r, m):
    """
       Function to compute sample entropy

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
       entropy: numpy.float64 (computed as -np.log(A / B))
           approximate entropy

       Example
       -------
       sample_entropy = approximate_entropy(data, 0.1*np.nanstd(data))
    """
    N = sig.shape[0]

    xlen = N - m
    x = np.full((xlen, m), np.inf, dtype='float64')
    for i in range(N - m):
        x[i] = sig[i: i + m]

    x_B = np.full((xlen + 1, m), np.inf, dtype='float64')
    for i in range(xlen + 1):
        x_B[i] = sig[i: i + m]

    # Save all matches minus the self-match, compute B
    B = cnt = 0
    for x_i in x:
        for x_j in x_B:
            if _maxdist(x_i, x_j) <= r:
                cnt += 1
        B += cnt-1
        cnt = 0

    # Similar for computing A
    m += 1
    x_A = np.full((N - m + 1, m), np.inf, dtype='float64')
    for i in range(N - m + 1):
        x_A[i] = sig[i: i + m]

    A = cnt = 0
    for x_i in x_A:
        for x_j in x_A:
            if _maxdist(x_i, x_j) <= r:
                cnt += 1
        A += cnt - 1
        cnt = 0

    return -np.log(A / B)


class SampleEntropy(Method):

    algorithm = 'SAMPLE_ENTROPY'
    algorithm_type = 'univariate'
    version = '1.0.0'
    dtype = [('sampen', 'float32')]

    def __init__(self, **kwargs):
        """
        Sample entropy

        Parameters
        ----------
        sig: np.ndarray
            1D signal
        m: int
            window length of compared run of data, recommended (2-8)
        r: float64
            filtering treshold, recommended values: (0.1-0.25)*std
        """

        super().__init__(compute_sample_entropy, **kwargs)
        self._event_flag = False
