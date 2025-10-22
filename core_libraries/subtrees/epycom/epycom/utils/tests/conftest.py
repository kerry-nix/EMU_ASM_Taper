# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.

# Std imports

# Third pary imports
import pytest

import numpy as np

# Local imports


@pytest.fixture(scope="module")
def create_testing_data():
    """
    Creates testing data
    """

    freqs = [2.5, 6.0, 10.0, 16.0, 32.5, 67.5, 165.0, 425.0, 800.0, 1500.0]

    fs = 5000
    n = fs*10
    data = np.zeros(n)
    basic_amp = 10

    x = np.arange(n)
    for freq in freqs:
        freq_amp = basic_amp / freq
        y = np.sin(2 * np.pi * freq * x / fs)
        data += y

    return data
