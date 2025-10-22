# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.

# Std imports

# Third pary imports
import numpy as np

# Local imports
from epycom.validation.feature_evaluation import (eval_feature_differences,
                                                  get_feature_differences)

from epycom.validation.util import (match_detections,
                                    check_detection_overlap)

from epycom.validation.precision_recall import (calculate_f_score,
                                                create_precision_recall_curve)

# ----- Utils -----


def test_check_detection_overlap():
    detection_1 = [3, 8]

    # No overlap
    detection_2 = [1, 2]
    assert check_detection_overlap(detection_1, detection_2) == False

    # Left overlap
    detection_2 = [1, 4]
    assert check_detection_overlap(detection_1, detection_2) == True

    # Right overlap
    detection_2 = [7, 9]
    assert check_detection_overlap(detection_1, detection_2) == True

    # Detection 2 in detection 1
    detection_2 = [4, 6]
    assert check_detection_overlap(detection_1, detection_2) == True

    # Detection 1 in detection 2
    detection_2 = [2, 9]
    assert check_detection_overlap(detection_1, detection_2) == True


def test_match_detections(create_testing_dataframe):
    df_1 = create_testing_dataframe[0]
    df_2 = create_testing_dataframe[1]

    bn = ['event_start', 'event_stop']
    freq_name = 'frequency'

    match_df = match_detections(df_1, df_2, bn, freq_name,
                                sec_unit=None, sec_margin=1)

    assert np.all(match_df['gs_index'] == match_df['dd_index'])


# ----- Feature evaluation -----
def test_get_feature_differences(create_testing_dataframe):
    df_1 = create_testing_dataframe[0]
    df_2 = create_testing_dataframe[1]

    bn = ['event_start', 'event_stop']
    feature_names = {'frequency': 'frequency',
                     'amplitude': 'amplitude'}

    match_df, N_missed = get_feature_differences(df_1, df_2, bn, feature_names)

    # Check if indices are all the same
    assert np.all(match_df['gs_index'] == match_df['dd_index'])

    # Check that frequency differences are zero
    assert np.sum(match_df['frequency_diff'] == 0) == 2

    # Check that amplitude dfferences are zero
    assert np.sum(match_df['amplitude_diff'] == 0) == 2

    # Check that none of
    assert N_missed == 0


def test_eval_feature_differences(create_testing_dataframe):
    df_1 = create_testing_dataframe[0]
    df_2 = create_testing_dataframe[1]

    bn = ['event_start', 'event_stop']
    feature_names = {'frequency': 'frequency',
                     'amplitude': 'amplitude'}

    match_df, _ = get_feature_differences(df_1, df_2, bn, feature_names)

    feature_diff_keys = ['frequency_diff', 'amplitude_diff']

    diff_df = eval_feature_differences(match_df, feature_diff_keys)

    assert round(diff_df['frequency_diff'], 5) == 0.21517
    assert round(diff_df['amplitude_diff'], 5) == 0.21517


# ----- Precision - recall -----
def test_calculate_f_score(create_testing_dataframe):
    df_1 = create_testing_dataframe[0]
    df_2 = create_testing_dataframe[1]

    bn = ['event_start', 'event_stop']

    p, r, f1 = calculate_f_score(df_1, df_2, bn, sec_unit=None, sec_margin=1)

    assert p == 1.0
    assert r == 1.0
    assert f1 == 1.0


def test_create_precision_recall_curve(create_testing_dataframe):
    df_1 = create_testing_dataframe[0]
    df_2 = create_testing_dataframe[1]

    bn = ['event_start', 'event_stop']

    p, r, f1 = create_precision_recall_curve(df_1, df_2, bn, 'threshold',
                                             sec_unit=None, sec_margin=1,
                                             eval_type='equal')

    assert np.sum(p) == 2.0
    assert np.sum(r) == 1.0
    assert round(np.sum(f1), 5) == 1.33333
