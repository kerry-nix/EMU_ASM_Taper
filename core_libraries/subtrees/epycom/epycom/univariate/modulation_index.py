# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.

# Std imports

# Third pary imports
import numpy as np
import scipy.signal as sp

# Local imports
from ..utils.method import Method


def compute_mi_count(data, nbins=18):
    """
    Function to compute modulation index (MI) of given data

    Parameters
    ----------
    data: numpy.ndarray
        data from which MI is computed
    nbins: int
        number of bins in which data will be separated, can affecct the result, default is 18

    Returns
    -------
    MI: float64
        modulation index computed as KL/np.log(nbins)

    Example
    -------
    MI = compute_mi_count(data)

    """

    size = 2 * np.pi / nbins
    position = np.zeros(nbins)
    mean_amp = np.zeros(nbins)

    # Binning the phases
    for bins in range(0, nbins):
        position[bins] = -np.pi + bins * size

    f_data = sp.hilbert(data)
    ampl = np.abs(f_data)
    ph = np.angle(f_data)

    # Computing average amplitude
    for j in range(0, nbins):
        phases1 = ampl[np.where(position[j] <= ph)]
        phases2 = ampl[np.where(ph < position[j] + size)]
        phases = np.intersect1d(phases1, phases2)
        mean_amp[j] = np.mean(phases)

    # Normalizing amplitude
    p = mean_amp / np.sum(mean_amp)

    # Computing Shannon entropy
    H = -np.sum(p * np.log(p))

    # Computing Kullback–Leibler distance
    KL = np.log(nbins) - H

    return KL / np.log(nbins)


class ModulationIndex(Method):

    algorithm = 'MODULATION_INDEX'
    algorithm_type = 'univariate'
    version = '1.0.0'
    dtype = [('mi', 'float32')]

    def __init__(self, **kwargs):
        """
        Modulation Index

        Parameters
        ----------
        data: numpy.ndarray
            data from which MI is computed
        nbins: int
            number of bins in which data will be separated, can affecct the result, default is 18
        """

        super().__init__(compute_mi_count, **kwargs)
