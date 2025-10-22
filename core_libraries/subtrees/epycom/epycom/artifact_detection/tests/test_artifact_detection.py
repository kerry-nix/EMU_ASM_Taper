# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering;
# Institute of Scientific Instruments of the CAS, v. v. i., Medical signals -
# Computational neuroscience. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.

# Std imports
from math import isclose

# Third pary imports

# Local imports
from epycom.artifact_detection import Saturation, PowerlineNoise


def test_saturation(create_testing_artifact_data, benchmark):
    compute_instance = Saturation()
    res = benchmark(compute_instance.run_windowed,
                    create_testing_artifact_data, 50000)
    compute_instance.run_windowed(create_testing_artifact_data,
                                  5000,
                                  n_cores=2)

    expected_val = -8.895755358781173e-07

    assert isclose(res[0][0], expected_val, abs_tol=10e-8)
    
def test_powerline_noise(create_testing_artifact_data, benchmark):
    compute_instance = PowerlineNoise(fs=5000)
    res = benchmark(compute_instance.run_windowed,
                     create_testing_artifact_data, 50000)
    compute_instance.run_windowed(create_testing_artifact_data,
                                  5000,
                                  n_cores=2)

    expected_val = 0.006305294

    assert isclose(res[0][0], expected_val, abs_tol=10e-5)
