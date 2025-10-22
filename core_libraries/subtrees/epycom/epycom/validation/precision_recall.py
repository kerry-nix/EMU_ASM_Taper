# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.


# Std imports

# Third pary imports

# Local imports
from epycom.validation.util import check_detection_overlap

"""
NOTE: we could use scikit-learn for this but that would require additional
modules to be installed. Something to consider for the future. Would simplify
things a bit. If we incorporate clustering and machine learning we should
switch to this.
"""


def create_precision_recall_curve(gs_df, dd_df, bn, threshold,
                                  sec_unit=None, sec_margin=1,
                                  eval_type='equal'):
    """
    Function to create precision recall curve.

    Parameters
    ----------
    gs_df: pandas.DataFrame
        Gold standard detections
    dd_df: pandas.DataFrame
        Automatically detected detections
    bn: list
        Names of event start stop [start_name, stop_name]
    threshold: str
        Name of the threshold field for evaluation
    sec_unit: int
        Nnumber representing one second of signal - this can
        significantly imporove the speed of this operation
    sec_margin: int
        Margin for creating subsets of compared data - should be set according
        to the legnth of compared events (1s for HFO should be enough)
    eval_type: str
        Whether to use bigger than threshold or equal to threshold
        for thresholding, options are 'equal' or 'bigger'

    Returns
    -------
    precision: list
        List of precision points
    recall: list
        List of recall points
    f1_score:
        List of f1 points
    """

    # Initiate lists
    precision = []
    recall = []
    f1_score = []

    # Thresholds
    ths = list(dd_df[threshold].unique())
    ths.sort()

    # Run through thresholds
    for th in ths:
        print('Processing threshold ' + str(th))

        if eval_type == 'equal':
            sub_dd_df = dd_df[dd_df[threshold] == th].copy()
        elif eval_type == 'bigger':
            sub_dd_df = dd_df[dd_df[threshold] >= th].copy()
        else:
            raise RuntimeError('Unknown eval_type "' + eval_type + '"')

        if sec_unit:
            p, r, f = calculate_f_score(gs_df,
                                        sub_dd_df,
                                        bn,
                                        sec_unit,
                                        sec_margin)
        else:
            p, r, f = calculate_f_score(gs_df,
                                        sub_dd_df,
                                        bn)

        precision.append(p)
        recall.append(r)
        f1_score.append(f)

    return precision, recall, f1_score


def calculate_f_score(gs_df, dd_df, bn, sec_unit=None, sec_margin=1):
    """
    Function to calculate precision and recall values.

    Parameters
    ----------
    gs_df: pandas.DataFrame
        Gold standard detections
    dd_df: pandas.DataFrame
        Automatically detected detections
    bn: list
        Names of event start stop [start_name, stop_name]
    threshold: str
        Name of the threshold field for evaluation
    sec_unit: int
        Nnumber representing one second of signal - this can
        significantly imporove the speed of this operation
    sec_margin: int
        Margin for creating subsets of compared data - should be set according
        to the legnth of compared events (1s for HFO should be enough)

    Returns
    -------
    precision: list
        List of precision points
    recall: list
        List of recall points
    f1_score:
        List of f1 points
    """

    # Create column for matching
    dd_df.loc[:, 'match'] = False
    gs_df.loc[:, 'match'] = False

    # Initiate true positive
    TP = 0

    # Start running through gold standards
    for gs_row in gs_df.iterrows():
        gs_det = [gs_row[1][bn[0]], gs_row[1][bn[1]]]

        det_flag = False
        if sec_unit:
            subset = dd_df[(dd_df[bn[0]] < gs_det[0] + sec_unit * sec_margin) &
                           (dd_df[bn[0]] > gs_det[0] - sec_unit * sec_margin)]
            for dd_row in subset.iterrows():
                dd_det = [dd_row[1][bn[0]], dd_row[1][bn[1]]]

                if check_detection_overlap(gs_det, dd_det):
                    det_flag = True
                    break
        else:
            for dd_row in dd_df.iterrows():
                dd_det = [dd_row[1][bn[0]], dd_row[1][bn[1]]]

                if check_detection_overlap(gs_det, dd_det):
                    det_flag = True
                    break

        # Mark the detections
        if det_flag:
            TP += 1
            dd_df.loc[dd_row[0], 'match'] = True
            gs_df.loc[gs_row[0], 'match'] = True

    # We ge number of unmatched detections
    FN = len(gs_df[gs_df['match'] == False])
    FP = len(dd_df[dd_df['match'] == False])

    # Calculate precision and recall
    precision = TP / (TP + FP)
    recall = TP / (TP + FN)
    if precision == 0 and recall == 0:
        f1_score = 0
    else:
        f1_score = 2 * ((precision * recall) / (precision + recall))

    return precision, recall, f1_score
