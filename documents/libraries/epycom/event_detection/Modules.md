# epycom.event_detection.BarkmeierDetector


=======

# epycom.event_detection.CSDetector


=======

# epycom.event_detection.HilbertDetector


=======

# epycom.event_detection.LineLengthDetector


=======

# epycom.event_detection.RootMeanSquareDetector


=======

# epycom.event_detection.detect_hfo_cs_beta

    Beta version of CS detection algorithm. Which was used to develop
    CS detection algorithm.

    Parameters
    ----------
    sig: numpy array
        1D numpy array with raw data
    fs: int
        Signal sampling frequency
    threshold: float
        Threshold for detection between 0 and 1 (Default=0.1)
    cycs_per_detect: float
        Minimal number of cycles threshold. (Default=4)
    mp: int
        Number of cores to use (def = 1)

    Returns
    -------
    output: list
        List of tuples with the following structure of detections:
        (event_start, event_stop, low_fc, high_fc, amp, fhom, dur, prod,
         is_conglom)

    Note
    ----
    The last object "is_conglom" in the output bool value.
    False stands for detections in single frequency bands whereas True
    stands for conglomerate created from detections in frequency bands that
    overlap in time domain.

    References
    ----------
    [1] J. Cimbalnik, A. Hewitt, G. A. Worrell, and M. Stead, “The CS
    Algorithm: A Novel Method for High Frequency Oscillation Detection
    in EEG,” J. Neurosci. Methods, vol. 293, pp. 6–16, 2017.
    

=======

# epycom.event_detection.detect_hfo_hilbert

    Slightly modified algorithm which uses the 2D HFO hilbert detection
    used in Kucewicz et al. 2014.

    Parameters
    ----------
    sig: numpy array
        1D numpy array with EEG data
    fs: float
        Sampling frequency of the signal
    low_fc: float
        Low cut-off frequency
    high_fc: float
        High cut-off frequency
    threshold: float
        Threshold for detection (default=3)
    band_spacing: str
        Spacing of hilbert freqeuncy bands - options: 'linear' or 'log'
        (default='linear'). Linear provides better frequency resolution but
        is slower.
    num_bands: int
        Number of bands if band_spacing = log (default=300)
    cyc_th: float
        Minimum number of cycles to detect (deafult=1)
    gap_th: float
        Number of cycles for gaps (default=1)
    mp: int
        Number of cores to use (default=1)

    Returns
    -------
    output: list
        List of tuples with the follwoing structure of detections:
        (event_start, event_stop, freq_min, freq_max, freq_at_max,
         max_amplitude)

    References
    ----------
    [1] M. T. Kucewicz, J. Cimbalnik, J. Y. Matsumoto, B. H. Brinkmann,
    M. Bower, V. Vasoli, V. Sulc, F. Meyer, W. R. Marsh, S. M. Stead, and
    G. A. Worrell, “High frequency oscillations are associated with cognitive
    processing in human recognition memory.,” Brain, pp. 1–14, Jun. 2014.
    

=======

# epycom.event_detection.detect_hfo_ll

    Line-length detection algorithm.

    Parameters
    ----------
    sig: np.ndarray
        1D array with raw data (already filtered if required)
    fs: int
        Sampling frequency
    threshold: float
        Number of standard deviations to use as a threshold
    window_size: int
        Sliding window size in samples
    window_overlap: float
        Fraction of the window overlap (0 to 1)

    Returns
    -------
    output: list
        List of tuples with the following structure:
        (event_start, event_stop)

    References
    ----------
    [1] A. B. Gardner, G. A. Worrell, E. Marsh, D. Dlugos, and B. Litt, “Human
    and automated detection of high-frequency oscillations in clinical
    intracranial EEG recordings,” Clin. Neurophysiol., vol. 118, no. 5, pp.
    1134–1143, May 2007.
    

=======

# epycom.event_detection.detect_hfo_rms

    Root mean square detection algorithm {Staba et al. 2002,
    Blanco et al 2010}.

    Parameters
    ----------
    sig: np.ndarray
        1D array with raw data (already filtered if required)
    fs: int
        Sampling frequency
    threshold: float
        Number of standard deviations to use as a threshold
    window_size: int
        Sliding window size in samples
    window_overlap: float
        Fraction of the window overlap (0 to 1)

    Returns
    -------
    output: list
        List of tuples with the following structure:
        (event_start, event_stop)

    References
    ----------
    [1] R. J. Staba, C. L. Wilson, A. Bragin, I. Fried, and J. Engel,
    “Quantitative Analysis of High-Frequency Oscillations (80 − 500 Hz)
    Recorded in Human Epileptic Hippocampus and Entorhinal Cortex,”
    J. Neurophysiol., vol. 88, pp. 1743–1752, 2002.
    

=======

# epycom.event_detection.detect_spikes_barkmeier

    Python version of Barkmeier's EEG spike detector. {Barkmeier et al. 2011}

    Parameters
    ----------
    sig: np.ndarray
        1D numpy array of EEG data
    fs: int
        sampling frequency of the signal
    scale: float        scaling parameter (def=70)
    std_coeff: float
        z-score threshold for spike detection (def=4)
    through_search: float
        extent to which search for spike throughs in s (def=0.04)
    det_thresholds: dict
        detection thresholds (dictionary)
        {'LS':700, # Left slope
         'RS':700, # Right slope
         'TAMP':600, # Total amplitude
         'LD':0.01, # Left duration
         'RD':0.01} # Right duration
    filter_spec: dict
        narrow and broad band filter specifications
        {'narrow':[20, 50],
         'broad':[1, 80]}
    win_idx: int
        Statistical window index. This is used when the
        function is run in separate windows. Default = None

    Returns
    -------
    output: list
        List of tuples with the following structure of detections:
        (event_peak, event_amp, left_amp, left_dur, right_amp, right_dur)
    

=======

# epycom.event_detection.hfo


=======

# epycom.event_detection.spike


=======

