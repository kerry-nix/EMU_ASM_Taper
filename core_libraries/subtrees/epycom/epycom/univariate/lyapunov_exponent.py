# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.

# Std imports

# Third pary imports
import numpy as np
from scipy.spatial.distance import pdist, squareform

# Local imports
from ..utils.method import Method


def _compute_phase_space(data, dimensions, sample_lag):
    """
    Create phase space of time series data. This function takes value lagged by
     sample_lag as coordination in new dimension.

    Parameters
    ----------
    data: np.array
        Signal to analyze, time series (array, int, float)
    dimensions: int
        Number of dimensions to create (int)
    sample_lag: int
        Delay in samples used for coordination extraction (int)

    Returns
    -------
    space: numpy.ndarray
        "Dimensions" times dinmensional space

    Example
    -------
    space = _compute_phase_space(data,8,5000)
    """

    length = np.size(data)
    space = np.zeros([dimensions, length - ((dimensions - 1) * sample_lag)])
    for i in range(dimensions):
        start = i * sample_lag
        end = length - ((dimensions - i - 1) * sample_lag)
        space[i, :] = data[start:end]
    return space


def _compute_acorr_exp(data, fs):
    """
    Find point, where autocorrelation drops to 1-1/np.e of it's maximum

    Paremeters
    ----------
    data: np.array
        Signal to analyze, time series (array, int, float)
    fs: float
        Sampling frequency

    Returns
    -------
    point: sample where autocorrelation drops to  1-1/np.e of it's
            maximum
    """

    # It is supposed that autocorrelation function of EEG drops to 1-1/e of
    # it's value within one second. This assumption masively reduces
    # computation time
    data = data[0:fs]

    # normalize data
    data = data - np.mean(data)

    acorr = np.correlate(data, data, mode='full')
    acorr = acorr / max(acorr)
    acorr = acorr[acorr > (1 - 1 / np.e)]
    point = len(acorr) // 2

    return point


def compute_lyapunov_exponent(data, fs=5000, dimension=5, sample_lag=None,
                              trajectory_len=20, min_tsep=500):
    """
    Lyapnov largest exponent estimation according to Rosenstein algorythm

    With use of some parts from nolds library:
    https://pypi.org/project/nolds
    https://github.com/CSchoel

    Parameters
    ----------
    data: np.array
        Signal to analyze, time series (array, int, float).
    fs: float
        Sampling frequency
    dimension: int
        Number of dimensions to compute lyapunov exponent.
    sample_lag: int
        Delay in samples used for coordination extraction.
    trajectory_len: int
        Number of points on divergence trajectory.
    min_tsep: int
        Nearest neighbors have temporal separation greater then min_tstep.

    Returns
    -------
    le: float
        Estimation of largest Lyapunov coeficient acording to Rosenstein
        algorithm.

    Example
    -------
    le = compute_lyapunov_exp(data, fs=5000, dimension=5, sample_lag=None,
                         trajectory_len=20, min_tsep=500)
    """

    # If sample lag for creating orbit is not set, it will be counted as
    # a point, where autocorrelation function is 1-1/e.
    if sample_lag is None:
        sample_lag = _compute_acorr_exp(data, fs)

    # creating m-dimensional orbit by delaying 1D signal
    orbit = _compute_phase_space(data, int(dimension), int(sample_lag))

    # calculate euclidian distances amog all points in orbit
    distances = squareform(pdist(orbit.T, 'euclidean'))

    m = len(distances)

    # we do not want to consider vectors as neighbor that are less than
    # min_tsep time steps together => mask the distances min_tsep to the right
    # and left of each index by setting them to infinity (will never be
    # considered as nearest neighbors)

    for i in range(m):
        distances[i, max(0, i - min_tsep):i + min_tsep + 1] = float("inf")

    # check that we have enough data points to continue
    ntraj = m - trajectory_len + 1
    min_traj = min_tsep * 2 + 2  # in each row min_tsep + 1 disances are inf
    if ntraj <= 0:
        msg = "Not enough data points. Need {} additional data points to " \
            + "follow a complete trajectory."
        raise ValueError(msg.format(-ntraj + 1))
    if ntraj < min_traj:
        # not enough data points => there are rows where all values are inf
        assert np.any(np.all(np.isinf(distances[:ntraj, :ntraj]), axis=1))
        msg = "Not enough data points. At least {} trajectories are " \
            + "required to find a valid neighbor for each orbit vector with " \
            + "min_tsep={} but only {} could be created."
        raise ValueError(msg.format(min_traj, min_tsep, ntraj))
    assert np.all(np.any(np.isfinite(distances[:ntraj, :ntraj]), axis=1))

    # find nearest neighbors (exclude last columns, because these vectors
    # cannot be followed in time for trajectory_len steps)
    nb_idx = np.argmin(distances[:ntraj, :ntraj], axis=1)

    # build divergence trajectory by averaging distances along the trajectory
    # over all neighbor pairs
    div_traj = np.zeros(trajectory_len, dtype=float)
    for k in range(trajectory_len):
        # calculate mean trajectory distance at step k
        indices = (np.arange(ntraj) + k, nb_idx + k)
        div_traj_k = distances[indices]
        # filter entries where distance is zero (would lead to -inf after log)
        nonzero = np.where(div_traj_k != 0)
        if len(nonzero[0]) == 0:
            # if all entries where zero, we have to use -inf
            div_traj[k] = -np.inf
        else:
            div_traj[k] = np.mean(np.log(div_traj_k[nonzero]))

    # filter -inf entries from mean trajectory
    ks = np.arange(trajectory_len)
    finite = np.where(np.isfinite(div_traj))
    ks = ks[finite]
    div_traj = div_traj[finite]

    if len(ks) < 1:
        # if all points or all but one point in the trajectory is -inf, we
        # cannot fit a line through the remaining points => return -inf as
        # exponent
        poly = [-np.inf, 0]
    else:
        # normal line fitting
        poly = np.polyfit(np.arange(len(div_traj)), div_traj, 1)

    le = poly[0] / (sample_lag / fs)

    return le


class LyapunovExponent(Method):

    algorithm = 'LYAPUNOV_EXPONENT'
    algorithm_type = 'univariate'
    version = '1.0.0'
    dtype = [('lyapunov_exponent', 'float32')]

    def __init__(self, **kwargs):
        """
        Lyapnov largest exponent estimation according to Rosenstein algorythm

        With use of some parts from nolds library:
        https://pypi.org/project/nolds
        https://github.com/CSchoel

        Parameters
        ----------
        fs: float
            Sampling frequency
        dimensions: int
            Number of dimensions to compute lyapunov exponent.
        sample_lag: int
            Delay in samples used for coordination extraction.
        trajectory_len: int
            Number of points on divergence trajectory.
        min_tstep: int
            Nearest neighbors have temporal separation greater then min_tstep.
        """

        super().__init__(compute_lyapunov_exponent, **kwargs)
