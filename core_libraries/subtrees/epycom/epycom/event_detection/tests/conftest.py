# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.

# Std imports

# Third pary imports
import pytest

import numpy as np

# Local imports
from epycom.simulation.create_simulated import simulate_hfo, simulate_spike


@pytest.fixture(scope="module")
def create_testing_eeg_data():
    """
    Creates testing data
    """

    freqs = [2.5, 6.0, 10.0, 16.0, 32.5, 67.5, 165.0,
             250.0, 425.0, 500.0, 800.0, 1500.0]

    fs = 5000
    n = fs * 10
    data = np.zeros(n)
    basic_amp = 10

    x = np.arange(n)
    for freq in freqs:
        freq_amp = basic_amp / freq
        y = np.sin(2 * np.pi * freq * x / fs)
        data += y

    # We have dummy data now inject 2 HFOs and a spike
    fs = 5000
    freq = 250
    numcycles = 9
    sim = simulate_hfo(fs, freq, numcycles)[0]
    ev_start = 5000
    data[ev_start: ev_start+len(sim)] += sim*10

    fs = 5000
    dur = 0.1
    sim = simulate_spike(fs, dur)
    ev_start = 4*5000
    data[ev_start: ev_start+len(sim)] += sim*30

    fs = 5000
    freq = 500
    numcycles = 9
    sim = simulate_hfo(fs, freq, numcycles)[0]
    ev_start = 7*5000
    data[ev_start: ev_start+len(sim)] += sim*10

    return data
