# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.


# Std imports

# Third pary imports
from scipy.stats import ttest_1samp

# Local imports
from epycom.validation.util import match_detections


def eval_feature_differences(diff_df, feature_diff_keys):
    """
    Function to evaluate feature differences between known values and estimated
    values.

    Parameters
    ----------
    diff_df: pandas.DataFrame
        Dataframe produced by get_feature_differences
    feature_diff_keys: list
        Names columns to evaluate

    Returns
    -------
    stat_dict: dict
        Dictionary with p values of 1 sample ttests for each feature
    """

    # Run statistical test on
    stat_dict = {}
    for f_key in feature_diff_keys:
        res = ttest_1samp(diff_df.loc[:, f_key].values, 0)[1]
        stat_dict[f_key] = res

    return stat_dict


def get_feature_differences(gs_df, dd_df, bn, feature_names):
    """
    Function to get feature differences between known and estimated values.

    Parameters
    ----------
    gs_df: pandas.DataFrame
        Dataframe of events with known features
    dd_df: pandas.DataFrame
        Dataframe of events with estimated efatures
    bn: list
        Names of event start stop [start_name, stop_name]
    feature_names: dict
        Dictionary with features as keys and column names as values

    Returns
    -------
    match_df: pandas.DataFrame
        Dataframe with indexes of matched detections
    N_missed: int
        Number of missed gold standard detections
    """

    # Match the detections first
    if 'frequency' in feature_names.keys():
        match_df = match_detections(
            gs_df, dd_df, bn, feature_names['frequency'])
    else:
        match_df = match_detections(gs_df, dd_df, bn)

    # Get count of missed detections and pop them from the df
    N_missed = len(match_df.loc[match_df.dd_index.isnull()])
    match_df = match_df.loc[~match_df.dd_index.isnull()]

    # Run through matched detections and compare the features
    for feature in feature_names.keys():
        for match_row in match_df.iterrows():
            gs_feat = gs_df.loc[match_row[1].gs_index, feature_names[feature]]
            dd_feat = dd_df.loc[match_row[1].dd_index, feature_names[feature]]
            match_df.loc[match_row[0], feature + '_diff'] = gs_feat - dd_feat

    return match_df, N_missed
