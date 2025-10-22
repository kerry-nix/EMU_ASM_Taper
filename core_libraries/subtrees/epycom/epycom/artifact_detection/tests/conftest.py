# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering;
# Institute of Scientific Instruments of the CAS, v. v. i., Medical signals -
# Computational neuroscience. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.

# Std imports

# Third pary imports
import pytest

import numpy as np

# Local imports


@pytest.fixture(scope="module")
def create_testing_artifact_data():
    """
    Creates testing data
    """

    freqs = [2.5, 6.0, 10.0, 16.0, 32.5, 67.5, 165.0, 425.0, 800.0, 1500.0]
    amps = [1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.01, 0.001, 0.0005, 0.0001]
    
    fs = 5000
    n = fs*10
    data = np.zeros(n)

    x = np.arange(n)
    for i,freq in enumerate(freqs):
        a = amps[i]
        y = a * np.sin(2 * np.pi * freq * x / fs)
        data += y
                
    data[int(n/10):int(2*n/10)] = 0
    data[int(3*n/10):int(4*n/10)] = max(data)
    data[int(5*n/10)] = max(data)
    
    hf_noise = np.zeros(int(n/10))
    x = np.arange(int(n/10))
    amps = [0.3, 0.1, 0.01]
    for i,freq in enumerate([650.0, 800.0, 950.0]):
        a = amps[i]
        y = a * np.sin(2 * np.pi * freq * x / fs)
        hf_noise += y
    data[int(6*n/10):int(7*n/10)] = data[int(6*n/10):int(7*n/10)]+hf_noise
    
    x = np.arange(int(n/10))
    line_noise = 1 * np.sin(2 * np.pi * 50.0 * x / fs)
    data[int(8*n/10):int(9*n/10)] = data[int(8*n/10):int(9*n/10)]+line_noise

    return data
