# CNT_research_tools.bandpower
Adapted from https://raphaelvallat.com/bandpower.html
    Compute the average power of the signal x in a specific frequency band.

    Parameters
    ----------
    data : 1d-array or 2d-array
        Input signal in the time-domain. (time by channels)
    sf : float
        Sampling frequency of the data.
    band : list
        Lower and upper frequencies of the band of interest.
    window_sec : float
        Length of each window in seconds.
        If None, window_sec = (1 / min(band)) * 2
    relative : boolean
        If True, return the relative power (= divided by the total power of the signal).
        If False (default), return the absolute power.

    Return
    ------
    bp : float
        Absolute or relative band power.
    

=======

# CNT_research_tools.clean_labels

    This function cleans a list of channels and returns the new channels
    

=======

# CNT_research_tools.find_non_ieeg

    This function finds non-iEEG channel labels
    

=======

# CNT_research_tools.get_iEEG_data
"
    2020.04.06. Python 3.7
    Andy Revell, adapted by Akash Pattnaik (2021.06.23)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Purpose:
    To get iEEG data from iEEG.org. Note, you must download iEEG python package from GitHub - instructions are below
    1. Gets time series data and sampling frequency information. Specified electrodes are removed.
    2. Saves as a pickle format
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Input
        username: your iEEG.org username
        password_bin_file: your iEEG.org password bin_file
        iEEG_filename: The file name on iEEG.org you want to download from
        start_time_usec: the start time in the iEEG_filename. In microseconds
        stop_time_usec: the stop time in the iEEG_filename. In microseconds.
            iEEG.org needs a duration input: this is calculated by stop_time_usec - start_time_usec
        ignore_electrodes: the electrode/channel names you want to exclude. EXACT MATCH on iEEG.org. Caution: some may be LA08 or LA8
        outputfile: the path and filename you want to save.
            PLEASE INCLUDE EXTENSION .pickle.
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Output:
        Saves file outputfile as a pickle. For more info on pickling, see https://docs.python.org/3/library/pickle.html
        Briefly: it is a way to save + compress data. it is useful for saving lists, as in a list of time series data and sampling frequency together along with channel names
        List index 0: Pandas dataframe. T x C (rows x columns). T is time. C is channels.
        List index 1: float. Sampling frequency. Single number
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Example usage:
    username = 'arevell'
    password = 'password'
    iEEG_filename='HUP138_phaseII'
    start_time_usec = 248432340000
    stop_time_usec = 248525740000
    removed_channels = ['EKG1', 'EKG2', 'CZ', 'C3', 'C4', 'F3', 'F7', 'FZ', 'F4', 'F8', 'LF04', 'RC03', 'RE07', 'RC05', 'RF01', 'RF03', 'RB07', 'RG03', 'RF11', 'RF12']
    outputfile = '/Users/andyrevell/mount/DATA/Human_Data/BIDS_processed/sub-RID0278/eeg/sub-RID0278_HUP138_phaseII_248432340000_248525740000_EEG.pickle'
    get_iEEG_data(username, password, iEEG_filename, start_time_usec, stop_time_usec, removed_channels, outputfile)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    To run from command line:
    python3.6 -c 'import get_iEEG_data; get_iEEG_data.get_iEEG_data("arevell", "password", "HUP138_phaseII", 248432340000, 248525740000, ["EKG1", "EKG2", "CZ", "C3", "C4", "F3", "F7", "FZ", "F4", "F8", "LF04", "RC03", "RE07", "RC05", "RF01", "RF03", "RB07", "RG03", "RF11", "RF12"], "/gdrive/public/DATA/Human_Data/BIDS_processed/sub-RID0278/eeg/sub-RID0278_HUP138_phaseII_D01_248432340000_248525740000_EEG.pickle")'
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #How to get back pickled files
    with open(outputfile, 'rb') as f: data, fs = pickle.load(f)
    

=======

# CNT_research_tools.gini
Calculate the Gini coefficient of a numpy array.
    From https://github.com/oliviaguest/gini 
    

=======

# CNT_research_tools.line_length


=======

# CNT_research_tools.movmean


=======

# CNT_research_tools.plot_iEEG_data
"
    2021.06.23. Python 3.8
    Akash Pattnaik
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Purpose:
    To plot iEEG data
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Input
        data: iEEG data in pandas.DataFrame or numpy.array
        time: time array 
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Output:
        Returns figure handle
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    

=======

# CNT_research_tools.pull_patient_localization


=======

# CNT_research_tools.pull_sz_ends


=======

# CNT_research_tools.pull_sz_starts


=======

