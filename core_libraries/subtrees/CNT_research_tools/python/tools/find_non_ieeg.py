'''
    This function finds non-iEEG channel labels
'''
import re
import numpy as np

non_ieeg = [
    "EKG",
    "O",
    "C",
    "ECG"
    ]

def find_non_ieeg(channel_li):
    '''
    This function finds non-iEEG channel labels
    '''

    is_non_ieeg = np.zeros(len(channel_li), dtype=bool)
    for ind, i in enumerate(channel_li):
        # Gets letter part of channel
        regex_match = re.match(r"(\D+)(\d+)", i)
        lead = regex_match.group(1)

        # finds non-iEEG channels, make a separate function
        if lead in non_ieeg:
            is_non_ieeg[ind] = True

    return is_non_ieeg
    