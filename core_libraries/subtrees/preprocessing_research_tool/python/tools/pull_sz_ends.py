import numpy as np


def pull_sz_ends(patient, metadata):
    assert patient in metadata
    sz_names = metadata[patient]["Events"]["Ictal"]

    sz_ends = []
    for i_sz, sz_name in enumerate(sz_names):
        if patient == "HUP111":
            if "D01" in sz_names[sz_name]["iEEG_record"]:
                continue
        if patient == "HUP181":
            if "D02" not in sz_names[sz_name]["iEEG_record"]:
                continue

        sz_ends.append(sz_names[sz_name]["SeizureEnd"])
    sz_ends = np.array(sz_ends)
    return np.unique(sz_ends)
