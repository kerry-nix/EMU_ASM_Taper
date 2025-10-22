import numpy as np
from scipy.spatial.distance import cdist

def laplacian_reference(data:np.ndarray,locs,radius,labels):
    """_summary_

    Args:
        data (np.ndarray): _description_
        locs (_type_): _description_
        radius (_type_): _description_
        labels (_type_): _description_
    """
    print("Currently only available for test for HUP212_phaseII")
    out_values = np.nan*np.zeros(data.shape)
    nchs = data.shape[1]

# Calculate pairwise distances between X and Y
    D = cdist(locs, locs)

    close = D < radius
    close[np.eye(close.shape[0]).astype(bool)] = 0

    close_chs = []

    for i in range(nchs):
        out_values[:,i] = data[:,i] - np.nanmean(data[:,close[i,:]],1)
        close_chs.append(np.where(close[i,:])[0])

    # laplacian_labels = [i+'-laplacian' for i in labels]
    laplacian_labels = []
    for i in range(len(labels)):
        ch_string = '-'
        ch_list = []
        for j in range(len(close_chs[i])):
            ch_list.append(labels[close_chs[i][j]])
        ch_string = ch_string + ",".join(ch_list)
        laplacian_labels.append(labels[i] + ch_string)

    return out_values,close_chs,laplacian_labels