#%%
import os
import sys

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from datetime import timedelta
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection

sys.path.insert(1, './../tools/')
from apple_parse import *
# %%
for subject in ['04','05','06','07','10','11']:
    lb3_id = "LB3_0" + subject + '_phaseII'
    path_to_lb3_folder = "/Users/jalpanchal/Documents/littlab/Data/"
    path_to_lb3_folder = "/gdrive/public/DATA/Human_Data/LB3_PIONEER/"
    path_to_subject_folder = os.path.join(path_to_lb3_folder,lb3_id)
    path_to_wearables = os.path.join(path_to_subject_folder,"wearables/")
    path_to_preprocessed = os.path.join(path_to_wearables, "pre-processed/")
    path_to_raw_activity = os.path.join(path_to_wearables,"raw/apple_health_export/")
    path_to_activity_csvs = os.path.join(path_to_preprocessed,"activity/")

    sleep = pd.read_csv(os.path.join(path_to_preprocessed, "watch_sleep.csv"))

    lines = []
    colors = []
    n_nights = len(sleep.creationDate.unique())

    for y, (creationDate, night) in enumerate(sleep.groupby('creationDate')):
        start_date = pd.to_datetime(night['timeStamp'].iloc[0]).floor('d')
        if pd.to_datetime(night['timeStamp'].iloc[0]).hour < 12:
            start_date = start_date - timedelta(1)

        for ind, sleep_segment in night.iterrows():
            x_vals = pd.date_range(sleep_segment['timeStamp'], sleep_segment['endDate'], freq="1s")
            x_vals = x_vals - start_date
            x_vals = x_vals.total_seconds() / (60 * 60)

            y_vals = np.ones(len(x_vals)) * y

            lines.append(np.column_stack((x_vals, y_vals)))
            colors.append(mcolors.to_rgba('tab:blue'))

        wake_times = {}
        wake_times['start'] = night['endDate'].iloc[:-1].values
        wake_times['end'] = night['timeStamp'][1:].values
        wake_times = pd.DataFrame(wake_times)

        for ind, wake_segment in wake_times.iterrows():
            # there is one patient that was lucky enough to be recorded literally during daylight savings time
            try:
                x_vals = pd.date_range(wake_segment['start'], wake_segment['end'], freq="1s")
            except:
                continue
            x_vals = x_vals - start_date
            x_vals = x_vals.total_seconds() / (60 * 60)

            y_vals = np.ones(len(x_vals)) * y

            lines.append(np.column_stack((x_vals, y_vals)))
            colors.append(mcolors.to_rgba('tab:orange'))

    fig, ax = plt.subplots()
    ax.set_xlim(20, 34)
    ax.set_ylim(-1, n_nights)
    nights = sleep.creationDate.unique()
    nights = nights - nights[0] + 1
    nights = [f"Night {i}" for i in nights]

    line_segments = LineCollection(
        lines,
        colors=colors,
        linewidths=(20),
        linestyle='solid'
        )
    ax.add_collection(line_segments)
    ax.invert_yaxis()
    xticks = ax.get_xticks()
    ax.set_xticklabels(["{:2d}:00".format(int(i)) for i in xticks % 24])
    ax.set_yticks(np.arange(n_nights), labels=nights)
    ax.set_title( "LB3_0" + subject)
    sns.despine()
    plt.show()
# %%
