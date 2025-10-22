# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.


# Std imports

# Third pary imports
import numpy as np

# Local imports
from epycom.utils.data_operations import (calculate_absolute_samples,
                                          create_output_df,
                                          add_metadata,
                                          process_matrix)
from epycom.utils.signal_transforms import (compute_hilbert_envelope,
                                            compute_hilbert_power,
                                            compute_teager_energy,
                                            compute_rms,
                                            compute_stenergy,
                                            compute_line_lenght,
                                            compute_stockwell_transform)
from epycom.utils.thresholds import th_std, th_tukey, th_percentile, th_quian


# ----- Data operations -----
def test_calculate_absolure_samples():
    test_arr = np.array([0, 1], dtype=[('win_idx', 'int32')])
    comp_arr = calculate_absolute_samples(test_arr, window_size=1000,
                                          overlap=0.5)
    assert np.all(np.array([0, 500]) == comp_arr)

def test_create_output_df():
    res_df = create_output_df(fields={'field_1': np.int32,
                                      'field_2': np.float32})
    expected_columns = ['event_start', 'event_stop', 'field_1', 'field_2']
    assert expected_columns == list(res_df.columns)

def test_add_metadata():
    res_df = create_output_df(fields={'field_1': np.int32,
                                      'field_2': np.float32})
    metadata = {'field_3': np.int32,
                'field_4': np.float32}
    add_metadata(res_df, metadata)
    expected_columns = ['event_start', 'event_stop',
                        'field_1', 'field_2', 'field_3', 'field_4']
    assert expected_columns == list(res_df.columns)
    
def test_process_matrix():
    matrix = np.zeros([5,5], dtype=[('test_val', 'float'),
                                ('win_idx', 'int')])
    matrix[:, 2]['test_val'] = 10
    matrix['win_idx'] = np.arange(5)
    
    count_sig = process_matrix(matrix, 'test_val', th_percentile, 95)
    
    assert np.all(np.array([0, 0, 1, 0, 0]) == count_sig)

# ----- Signal transforms -----
def test_compute_hilbert_envelope(create_testing_data):
    assert (round(np.sum(compute_hilbert_envelope(create_testing_data)), 5)
            == round(141021.90763537044, 5))


def test_compute_hilbert_power(create_testing_data):
    assert (round(np.sum(compute_hilbert_power(create_testing_data)), 5)
            == round(499812.84844509006, 5))


def test_compute_teager_energy(create_testing_data):
    assert (round(np.sum(compute_teager_energy(create_testing_data)), 5)
            == round(96410.92390890958, 5))


def test_compute_rms(create_testing_data):
    assert (round(np.sum(compute_rms(create_testing_data)), 5)
            == round(101737.24636480425, 5))


def test_compute_stenergy(create_testing_data):
    assert (round(np.sum(compute_stenergy(create_testing_data)), 5)
            == round(249993.5292416787, 5))


def test_compute_line_lenght(create_testing_data):
    assert (round(np.sum(compute_line_lenght(create_testing_data)), 5)
            == round(58084.721256107114, 5))


def test_compute_stockwell_transform(create_testing_data):
    s = compute_stockwell_transform(create_testing_data, 5000, 80, 600)[0]
    assert round(np.abs(np.sum(np.sum(s))), 5) == round(75000.00000000402, 5)


# ----- Thresholds -----
def test_th_std(create_testing_data):
    assert (round(th_std(create_testing_data, 3), 5)
            == round(6.708203932499344, 5))


def test_th_tukey(create_testing_data):
    assert (round(th_tukey(create_testing_data, 3), 5)
            == round(10.659619047273361, 5))


def test_th_percentile(create_testing_data):
    assert (round(th_percentile(create_testing_data, 75), 5)
            == round(1.5228027210391037, 5))


def test_th_quian(create_testing_data):
    assert (round(th_quian(create_testing_data, 3), 5)
            == round(6.777704219110832, 5))
