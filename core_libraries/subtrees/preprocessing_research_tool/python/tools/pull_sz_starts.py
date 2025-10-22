import numpy as np


def pull_sz_starts(patient, metadata):
    assert patient in metadata
    sz_names = metadata[patient]["Events"]["Ictal"]

    sz_starts = []
    for sz_name in sz_names:
        if patient == "HUP111":
            if "D01" in sz_names[sz_name]["iEEG_record"]:
                continue
        if patient == "HUP181":
            if "D02" not in sz_names[sz_name]["iEEG_record"]:
                continue

        sz_starts.append(sz_names[sz_name]["SeizureEEC"])
    sz_starts = np.array(sz_starts)
    return np.unique(sz_starts)
