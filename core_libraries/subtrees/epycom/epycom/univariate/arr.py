# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.

# Std imports

# Third pary imports
import numpy as np
import scipy
from math import nan
from numba import njit

# Local imports
from ..utils.method import Method


@njit("i8(i8)", cache=True)
def _next_regular(target):
    """
    Find the next regular number greater than or equal to target.
    Regular numbers are composites of the prime factors 2, 3, and 5.
    Also known as 5-smooth numbers or Hamming numbers, these are the optimal
    size for inputs to FFTPACK.

    Parameters
    ----------
    target: int
        Target number for finding a regular number (must be a positive integer).

    Returns
    --------
    match/target: int
            The next regular number greater than or equal to target.
    """
    if target <= 6:
        return target

    # Quickly check if it's already a power of 2
    if not (target & (target - 1)):
        return target

    match = np.inf  # Anything found will be smaller
    p5 = 1
    while p5 < target:
        p35 = p5
        while p35 < target:
            # Ceiling integer division, avoiding conversion to float
            # (quotient = ceil(target / p35))
            quotient = -(-target // p35)
            # Quickly find next power of 2 >= quotient
            exp = 1
            num = quotient - 1
            while num // 2 != 0:
                num = num // 2
                exp += 1
            p2 = 2 ** exp  # ((quotient - 1).bit_length())

            N = p2 * p35
            if N == target:
                return N
            elif N < match:
                match = N
            p35 *= 3
            if p35 == target:
                return p35
        if p35 < match:
            match = p35
        p5 *= 5
        if p5 == target:
            return p5
    if p5 < match:
        match = p5
    return match


def _acovf(x):
    """
    Estimate autocovariances. Recoded from statsmodels package.

    Parameters
    ----------
    x : array_like
        Time series data. Must be 1d.

    Returns
    -------
    acov: ndarray
        The estimated autocovariances.

    References
    -----------
    .. [1] Parzen, E., 1963. On spectral analysis with missing observations
           and amplitude modulation. Sankhya: The Indian Journal of
           Statistics, Series A, pp.383-392.
    """
    xo = x - x.mean()

    n = len(x)
    d = n * np.ones(2 * n - 1)

    acov = np.correlate(xo, xo, "full")[n - 1:] / d[n - 1:]
    return acov


def _corrmtx(x, m):
    """
    Correlation matrix

    This function is used by PSD estimator functions. It generates
    the correlation matrix from a correlation data set and a maximum lag.
    Recorded from spectrum package.

    Parameters:
    -----------
    x_input: array
        autocorrelation samples (1D)
    m: int
        the maximum lag (Depending on the choice of the method, the 
        correlation matrix has different sizes, but the number of rows is
        always m+1)

    Returns
    -------
    c:
        the autocorrelation matrix

    Note
    ----
    Function implemented from original library spectrum
        https://pyspectrum.readthedocs.io/en/latest/#
    """

    from scipy.linalg import toeplitz
    N = len(x)

    return toeplitz(x[m:N], x[m::-1])


def _arcovar(x, order):
    """
    Simple and fast implementation of the covariance AR estimate.
    Recorded from spectrum package.

    Parameters
    ----------
    x: array
      Array of complex data samples
    order: int
        Order of linear prediction model

    Returns
    -------
    a: array
        Array of complex forward linear prediction coefficients
        (The output vector contains the normalized estimate of the AR system parameters)
    e: float64
        error

    Note
    ----
    Function implemented from original library spectrum
        https://pyspectrum.readthedocs.io/en/latest/#
    """

    X = _corrmtx(x, order)
    Xc = np.array(X[:, 1:])
    X1 = np.array(X[:, 0])

    # Coefficients estimated via the covariance method
    # Here we use lstsq rathre than solve function because Xc is not square
    # matrix

    a, _residues, _rank, _singular_values = scipy.linalg.lstsq(-Xc, X1)

    # Estimate the input white noise variance
    Cz = np.dot(X1.conj().transpose(), Xc)
    e = np.dot(X1.conj().transpose(), X1) + np.dot(Cz, a)
    assert e.imag < 1e-4, 'wierd behaviour'
    e = float(e.real)  # ignore imag part that should be small

    return a, e


@njit("f8(f8[:], i8)", cache=True)
def _arburg(x, order):
    """
      Estimate the complex autoregressive parameters by the Burg algorithm.

      .. math:: x(n) = qrt{(v}) e(n) + sum_{k=1}^{P+1} a(k) x(n-k)

      Parameters
      ----------
      x: numpy.ndarray
        array of complex data samples (length N)
      order: int
          order of autoregressive process (0<order<N)

      Returns
      -------
      rho: numpy.float64
          Real variable representing driving noise variance (mean square
          of residual noise) from the whitening operation of the Burg filter.

      Example
      -------
      v = arburg(sig, n)
      """
    N = len(x)

    # Initialisation
    # ------ rho, den
    rho = np.sum(np.abs(x) ** 2.) / float(N)  # Eq 8.21 [Marple]_
    den = rho * 2. * N

    a = np.zeros(order)

    ef = x.copy()
    eb = x.copy()

    temp = 1.

    #   Main recursion
    for k in range(0, order):
        num = 0
        # calculate the next order reflection coefficient Eq 8.14 Marple
        # numo = np.sum([ef[j] * eb[j - 1].conjugate() for j in range(k + 1, N)])
        for j in range(k + 1, N):
            num += ef[j] * eb[j - 1]

        den = temp * den - np.abs(ef[k]) ** 2 - np.abs(eb[N - 1]) ** 2
        kp = -2. * num / den  # eq 8.14

        temp = 1. - abs(kp) ** 2.
        new_rho = temp * rho

        # this should be after the criteria
        rho = new_rho
        if rho <= 0:
            return -1

        a_new = np.zeros(a.size + 1)
        for i in range(len(a)):
            a_new[i] = a[i]

        # a.resize(a.size + 1, refcheck=False)
        a_new[k] = kp
        if k == 0:
            for j in range(N - 1, k, -1):
                save2 = ef[j]
                ef[j] = save2 + kp * eb[j - 1]  # Eq. (8.7)
                eb[j] = eb[j - 1] + kp.conjugate() * save2

        else:
            # update the AR coeff
            khalf = int((k + 1) // 2)
            for j in range(0, khalf):
                ap = a_new[j]  # previous value
                a_new[j] = ap + kp * a_new[k - j - 1]  # Eq. (8.2)
                if j != k - j - 1:
                    a_new[k - j - 1] = a_new[k - j - 1] + kp * ap  # Eq. (8.2)

            # update the prediction error
            for j in range(N - 1, k, -1):
                save2 = ef[j]
                ef[j] = save2 + kp * eb[j - 1]  # Eq. (8.7)
                eb[j] = eb[j - 1] + kp.conjugate() * save2

    return rho


def _extract_poles(sig, n):
    """
    gives the complex poles from auto-regressive analysis

    Parameters
    ----------
    sig[channel, samples]: numpy.ndarray
    n: int
        requested order

    Returns
    -------
    roots: numpy.ndarray
        complex poles of the transfer function Z=exp(-u+2*Pii*w)
    v: numpy.float64
        residual variance(noise)
    A: numpy.ndarray
        complex amplitudes

    Example
    -------
     _, res, _ = _extract_poles(E[:, j], 3)
    """

    coeffs, _ = _arcovar(sig, n)  # complex forward linear prediction coefficients
    coeffs = np.append([1], coeffs)
    v = _arburg(sig, n)  # noise variance (mean square of residual noise)
    roots = np.roots(coeffs)
    R = np.transpose(roots * np.ones([roots.shape[0], roots.shape[0]]))
    R = R - np.transpose(R) + np.identity(roots.shape[0])
    A = 1 / np.prod(R, axis=0)

    return roots, v, A


def _embed(sig, tau, nED, step=1):
    """
    creates a time-delay embeded vector-signal from a 1D signal

    Parameters
    ----------
    sig: numpy.ndarray
        1D signal
    tau: int
        <-60> embedding delay, if tau<0 minimal correlation method is used in the range [1,-tau]
    nED: float
        number od embedding dimensions
    step: int
        default 1

    Returns
    -------
    V: list
        [embedding,time] vector signal
    tau: int
        the embedding delay (if automatically computed ), if <0 no zero-cross is detected

    Example
    -------
    E, _ = embed(sig, 1, 100.0, 50)
    """

    length = len(sig)
    tau1 = tau
    if tau < 0:
        tau = np.abs(tau)
        XC = _acovf(sig)  # estimated autocovariances
        XC = XC[1:]
        L = np.where((np.sign(XC[0:-1]) - np.sign(XC[1:])) != 0)  # indices of zero crossing
        if (len(L) > 0) and L[0] < tau:
            tau = L[0]
            tau1 = tau
        else:
            tau = np.argmin(np.abs(XC[0:tau])) + 1
            tau1 = -tau

    nv = np.floor((length - 1 - (nED - 1) * tau) / step) + 1
    V = [sig[0 + (j - 1) * tau::step][0:int(nv)] for j in range(1, int(nED + 1))]

    return V, tau1


def _get_residual(sig, n, winL):
    """
    Function calculates residuals

    Parameters
    ----------
    sig[channel,samples]: numpy.ndarray
    n: int
        requested model order
    winL: float
        window length with 50% overlap

    Returns
    -------
    res_var: list
        residual value
    """

    res_var = []

    E, _ = _embed(sig, 1, winL, round(winL / 2))
    E = np.array(E)

    for j in range(E.shape[1]):
        _, res, _ = _extract_poles(E[:, j], n)
        res_var.append(res)
        # roots.append(root)
        # ampt.append(A)

    # return np.array(roots).squeeze(), np.array(res_var), np.array(ampt).squeeze()
    return res_var


def compute_arr(sig, fs):
    """
    Function computes ARR parameters

    Parameters
    ----------
    sig: numpy.ndarray
    fs: float64
        sample frequency

    Returns
    -------
    ARRm: numpy.float64
        ARR parameters
    r1, r2, r3: list
                residuals for model order 1-3

    Example
    -------
    arrm = compute_arr(data, 5000)
    """

    # sig = stats.zscore(sig)
    winL = 0.02 * fs
    r1 = _get_residual(sig, 1, winL)
    r2 = _get_residual(sig, 2, winL)
    r3 = _get_residual(sig, 3, winL)

    # residual decline over residual from order 1 and 2 models
    rn = np.array([r1, r2])
    D = -2 * np.diff(rn, axis=0) / np.sum(rn, axis=0)
    r_clean = r3.copy()
    p95 = np.nanpercentile(r_clean, 95)
    w = D.shape[1] - 2

    while w > 0:
        if r_clean[w] > p95 and D[0][w] < 0.9:
            # remove possible artefacts at w and two neighbors windows
            # checking for border conditions
            r_clean[w] = nan  # is artefact

            if w != D.shape[1]:
                r_clean[w + 1] = nan

            if w != 0:
                r_clean[w - 1] = nan

            w = w - 1

        w = w - 1

    if any(sig[:]):
        # ARR = np.std(r3) / np.mean(r3)
        ARRm = np.nanstd(r_clean) / np.nanmean(r_clean)
    else:
        # ARR = nan
        ARRm = nan

    # return r1,r2,r3,r_clean,ARR,ARRm
    return ARRm


class AutoregressiveResidualModulation(Method):
    algorithm = 'AUTOROGRESSIVE_RESIDUAL_MODULATION'
    algorithm_type = 'univariate'
    version = '1.0.0'
    dtype = [('arr', 'float32')]

    def __init__(self, **kwargs):
        """
        Autoregressive residual modulation

        Parameters
        ----------
        sig:
            numpy.ndarray
        fsamp: float64
            sample frequency
        """

        super().__init__(compute_arr, **kwargs)
        self._event_flag = False
