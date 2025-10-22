import numpy as np
import pandas as pd

def automatic_bipolar_montage(data, data_columns):
    """This function returns the data in bipolar montage using the channel names

    Args:
        data (_type_): _description_
        data_columns (_type_): _description_

    Returns:
        _type_: _description_
    """
    channels = np.array(data_columns)

    nchan = len(channels)
    #naming to standard 4 character channel Name: (Letter)(Letter)[Letter](Number)(Number)
    # channels = channel2std(channels)
    count = 0
    for ch in range(nchan-1):
        ch1Ind = ch
        ch1 = channels[ch1Ind]
        #find sequential index
        ch2 = ch1[0:2] + f"{(int(ch1[2:4]) + 1):02d}"

        ch2exists = np.where(channels == ch2)[0]
        if len(ch2exists) > 0:
            ch2Ind = ch2exists[0]
            bipolar = pd.Series((data[:,ch1Ind] - data[:,ch2Ind])).rename(ch1)
            if count == 0: #initialize
                dfBipolar = pd.DataFrame( bipolar)
                count = count + 1
            else:
                dfBipolar = pd.concat([dfBipolar, pd.DataFrame( bipolar)], axis=1)
    return np.array(dfBipolar), np.array(dfBipolar.columns)