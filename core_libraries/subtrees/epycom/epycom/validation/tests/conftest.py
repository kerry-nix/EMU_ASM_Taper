# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.

# Std imports

# Third pary imports
import pytest

import pandas as pd
import numpy as np

# Local imports


@pytest.fixture(scope="module")
def create_testing_dataframe():
    """
    Creates testing data
    """

    data = np.array([[10, 20, 80, 5],
                     [30, 50, 100, 10],
                     [60, 70, 150, 15],
                     [80, 110, 200, 20]])

    dummy_df_1 = pd.DataFrame(data,
                              columns=['event_start', 'event_stop',
                                       'frequency', 'amplitude'])

    data = np.array([[10, 20, 80, 5, 1],
                     [30, 50, 100, 10, 1],
                     [60, 70, 152, 17, 2],
                     [80, 110, 201, 21, 2]])

    dummy_df_2 = pd.DataFrame(data,
                              columns=['event_start', 'event_stop',
                                       'frequency', 'amplitude', 'threshold'])

    

    return [dummy_df_1, dummy_df_2]