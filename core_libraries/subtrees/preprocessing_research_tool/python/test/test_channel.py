#%%
"""
This function provides a test case of pulling iEEG data
"""
# pylint: disable-msg=C0103
# pylint: disable-msg=C0301
#%%
# Imports
import pandas as pd
import numpy as np
import pytest
import os
import sys

work_path = os.getcwd()  # get current path
test_path = os.path.join(work_path, "python/test")
file_path = os.path.join(work_path, "python")
# print(os.getcwd())
# os.chdir(work_path)
# par_folder = os.path.dirname(work_path)
sys.path.append(file_path)

import tools

# %%
# unit test for clean_labels function
# write in a csv file, col 0 input, col 1 expected output
# wait to fetch all channel types from ieeg
test_chans = pd.read_csv(os.path.join(test_path, "decompLabel_testInput.csv"))
params = [
    tuple([[test_chans.iloc[i, 0]], [test_chans.iloc[i, 1]]])
    for i in range(test_chans.shape[0])
]


@pytest.mark.parametrize("input,expected", params)
def test_channel(input, expected):
    try:
        clean_channel = tools.clean_labels(input)
        assert expected == clean_channel
    except AttributeError as e:
        assert False
    # %%
