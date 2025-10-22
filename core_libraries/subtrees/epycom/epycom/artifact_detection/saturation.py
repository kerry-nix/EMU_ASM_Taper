# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering;
# Institute of Scientific Instruments of the CAS, v. v. i., Medical signals -
# Computational neuroscience. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.

# Std imports

# Third pary imports
import numpy as np

# Local imports
from ..utils.method import Method


def compute_saturation(sig):
    """
    Function to detect flat line (saturation or missing data)

    Parameters:
    ----------
    sig: np.array
        signal to analyze, time series (array, int, float)

    Returns
    -------
    results: tuple
        - mean_diff: average derivation in win

    Example
    -------
    sat = compute_saturation(sig)
    """

    return np.mean(np.diff(sig))


class Saturation(Method):

    algorithm = 'SATURATION'
    algorithm_type = 'artifact'
    version = '1.0.0'
    dtype = [('mean_diff', 'float32')]

    def __init__(self, **kwargs):
        """
        Class to detect flat line (saturation or missing data)

        Parameters:
        ----------
        sig: np.array
            signal to analyze, time series (array, int, float)

        """

        super().__init__(compute_saturation, **kwargs)
