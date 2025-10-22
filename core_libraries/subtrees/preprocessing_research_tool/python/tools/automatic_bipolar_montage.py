import numpy as np
import pandas as pd
import re
from beartype import beartype
from beartype.typing import Iterable

@beartype
def automatic_bipolar_montage(data: np.ndarray, labels: Iterable[str]):
    """This function returns the data in bipolar montage using the channel names

    Args:
        data (_type_): _description_
        labels (_type_): _description_

    Returns:
        _type_: _description_
    """
    channels = np.array(labels)

    nchan = len(channels)
    #chs_in_bipolar = np.empty(nchan,2)
    bipolar_labels = []
    out_values = data
    # naming to standard 4 character channel Name: (Letter)(Letter)[Letter](Number)(Number)
    # channels = channel2std(channels)
    for ch in range(nchan):
        out = np.nan*np.zeros(data.shape[0])
        ch1Ind = ch
        ch1 = channels[ch1Ind] # clean_label
        label_num_search = re.search(r"\d", ch1)
        if label_num_search is not None:
            label_num_idx = label_num_search.start()
            label_non_num = ch1[:label_num_idx]
            label_num = ch1[label_num_idx:]
            # find sequential index
            if int(label_num) > 12:
                print('This might be a grid and so bipolar might be tricky');
            ch2 = label_non_num + f"{int(label_num) + 1}"
            ch2exists = np.where(channels == ch2)[0]
            if len(ch2exists) > 0:
                ch2Ind = ch2exists[0]
                out = data[:,ch1Ind]-data[:,ch2Ind]
                bipolar_label = ch1+'-'+ch2
            else:
                bipolar_label = '-'
                #chs_in_bipolar(ch,:) = [ch1,ch2]
                #bipolar = pd.Series((data[:, ch1Ind] - data[:, ch2Ind])).rename(ch1)
                #if count == 0:  # initialize
                #    dfBipolar = pd.DataFrame(bipolar)
                #    count = count + 1
                #else:
                #    dfBipolar = pd.concat([dfBipolar, pd.DataFrame(bipolar)], axis=1)
        elif ch1 == 'FZ':
            ch2exists = np.where(channels == 'CZ')[0]
            if len(ch2exists) > 0:
                ch2Ind = ch2exists[0]
                out = data[:,ch1Ind]-data[:,ch2Ind]
                bipolar_label = ch1+'-'+'CZ'
                # chs_in_bipolar(ch,:) = [ch1,ch2]
        else:
            bipolar_label = '-'
        bipolar_labels.append(bipolar_label)
        out_values[:,ch] = out
    
    return out_values, np.array(bipolar_labels)

