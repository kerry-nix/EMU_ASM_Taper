# epycom.artifact_detection.PowerlineNoise


=======

# epycom.artifact_detection.Saturation


=======

# epycom.artifact_detection.compute_powerline_noise

    Function to detect the proportio of line noise in the signal

    Parameters:
    ----------
    sig: np.array
        signal to analyze, time series (array, int, float)
    fs: float
        sampling frequency
    freq: float or int
        line frequency (50 or 60 Hz)

    Returns
    -------
    score: float
        line frequency score

    Example
    -------
    sat = compute_saturation(sig)
    

=======

# epycom.artifact_detection.compute_saturation

    Function to detect flat line (saturation or missing data)

    Parameters:
    ----------
    sig: np.array
        signal to analyze, time series (array, int, float)

    Returns
    -------
    results: tuple
        - mean_diff: average derivation in win

    Example
    -------
    sat = compute_saturation(sig)
    

=======

# epycom.artifact_detection.powerline_noise


=======

# epycom.artifact_detection.saturation


=======

