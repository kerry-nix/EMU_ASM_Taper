# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.

"""
Some parts of code are recoded from package Anderson Brito da Silva's pyhfo
package (https://github.com/britodasilva/pyhfo)
"""

# Std imports

# Third pary imports
import numpy as np

# Local imports


def th_std(signal, ths):
    """
    Calcule threshold by Standar Desviations above the mean

    Parameters
    ----------
    signal: numpy array
        1D signal for threshold determination
    ths: float
        Number of SD above the mean

    Returns
    -------
    ths_value: float
        Value of the threshold
    """
    ths_value = np.mean(signal) + ths * np.std(signal)
    return ths_value


def th_tukey(signal, ths):
    """
    Calcule threshold by Tukey method.

    Parameters
    ----------
    signal: numpy array
        1D signal for threshold determination
    ths: float
        Number of interquartile interval above the 75th percentile

    Returns
    -------
    ths_value: float
        Value of the threshold
    """
    ths_value = np.percentile(
        signal, 75) + ths * (np.percentile(signal, 75)
                             - np.percentile(signal, 25))
    return ths_value


def th_percentile(signal, ths):
    """
    Calcule threshold by Percentile

    Parameters
    ----------
    signal: numpy array
        1D signal for threshold determination
    ths: float
        Percentile

    Returns
    -------
    ths_value: float
        Value of the threshold
    """
    ths_value = np.percentile(signal, ths)
    return ths_value


def th_quian(signal, ths):
    """
    Calcule threshold by Quian
    Quian Quiroga, R. 2004. Neural Computation 16: 1661–87.

    Parameters
    ----------
    signal: numpy array
        1D signal for threshold determination
    ths: float
        Number of estimated noise SD above the mean

    Returns
    -------
    ths_value: float
        Value of the threshold
    """
    ths_value = ths * np.median(np.abs(signal)) / 0.6745
    return ths_value
