# epycom.univariate.ApproximateEntropy


=======

# epycom.univariate.AutoregressiveResidualModulation


=======

# epycom.univariate.HjorthComplexity


=======

# epycom.univariate.HjorthMobility


=======

# epycom.univariate.LyapunovExponent


=======

# epycom.univariate.MeanVectorLength


=======

# epycom.univariate.ModulationIndex


=======

# epycom.univariate.PhaseLockingValue


=======

# epycom.univariate.PowerSpectralEntropy


=======

# epycom.univariate.SampleEntropy


=======

# epycom.univariate.ShannonEntropy


=======

# epycom.univariate.SignalStats


=======

# epycom.univariate.approximate_entropy


=======

# epycom.univariate.arr


=======

# epycom.univariate.compute_approximate_entropy

    Function computes approximate entropy of given signal

    Parameters
    ----------
    sig: np.ndarray
        1D signal
    r: np.float64
        filtering treshold, recommended values: (0.1-0.25)*np.nanstd(sig)
    m: int
        window length of compared run of data, recommended (2-8)

    Returns
    -------
    entro: numpy.float64
        approximate entropy

    Example
    -------
    signal_entropy = approximate_entropy(data, 0.1*np.nanstd(data))
    

=======

# epycom.univariate.compute_arr

    Function computes ARR parameters

    Parameters
    ----------
    sig: numpy.ndarray
    fs: float64
        sample frequency

    Returns
    -------
    ARRm: numpy.float64
        ARR parameters
    r1, r2, r3: list
                residuals for model order 1-3

    Example
    -------
    arrm = compute_arr(data, 5000)
    

=======

# epycom.univariate.compute_hjorth_complexity

    Compute Hjorth complexity of time series

    Parameters
    ----------
    signal: np.array
        Signal to analyze, time series (array, int, float)
    fs: float
        Sampling frequency of the time series

    Returns
    -------
    hjorth_complexity: float

    Example
    -------
    hjorth_complexity = compute_hjorth_complexity(data, 5000)

    Note
    ----
    result is NOT frequency dependent
    

=======

# epycom.univariate.compute_hjorth_mobility

    Function to compute Hjorth mobility of time series

    Parameters
    ----------
    signal: np.array
        Signal to analyze, time series (array, int, float)
    fs: float
        Sampling frequency of the time series

    Returns
    -------
    hjorth_mobility: float

    Example
    -------
    hjorth_mobility = compute_hjorth_mobility(data, 5000)

    Note
    ----
    Result is frequency dependent
    

=======

# epycom.univariate.compute_lyapunov_exponent

    Lyapnov largest exponent estimation according to Rosenstein algorythm

    With use of some parts from nolds library:
    https://pypi.org/project/nolds
    https://github.com/CSchoel

    Parameters
    ----------
    data: np.array
        Signal to analyze, time series (array, int, float).
    fs: float
        Sampling frequency
    dimension: int
        Number of dimensions to compute lyapunov exponent.
    sample_lag: int
        Delay in samples used for coordination extraction.
    trajectory_len: int
        Number of points on divergence trajectory.
    min_tsep: int
        Nearest neighbors have temporal separation greater then min_tstep.

    Returns
    -------
    le: float
        Estimation of largest Lyapunov coeficient acording to Rosenstein
        algorithm.

    Example
    -------
    le = compute_lyapunov_exp(data, fs=5000, dimension=5, sample_lag=None,
                         trajectory_len=20, min_tsep=500)
    

=======

# epycom.univariate.compute_mi_count

    Function to compute modulation index (MI) of given data

    Parameters
    ----------
    data: numpy.ndarray
        data from which MI is computed
    nbins: int
        number of bins in which data will be separated, can affecct the result, default is 18

    Returns
    -------
    MI: float64
        modulation index computed as KL/np.log(nbins)

    Example
    -------
    MI = compute_mi_count(data)

    

=======

# epycom.univariate.compute_mvl_count

    Function to compute mean vector lenght (MVL) of given data

    Parameters
    ----------
    fs: float64
        frequency
    data: numpy.ndarray
        data from which MI is computed
    lowband: list
            low frequency band boundaries [x, y], default [8, 12]
    highband: list
            high frequency band boundaries [x, y], default [250, 600]

    Returns
    -------
    mvl: numpy.complex128
        MVL of given signal

    Example
    -------
    MVL = compute_mvl_count(5000.0, data)

    

=======

# epycom.univariate.compute_plv_count

    Function to compute phase-locking value (PLV) of given data

    Parameters
    ----------
    fs: float64
        frequency
    data: numpy.ndarray
        data from which MI is computed
    lowband: list
            low frequency band boundaries [x, y], default [8, 12]
    highband: list
            high frequency band boundaries [x, y], default [250, 600]

    Returns
    -------
    MI: numpy.complex128
        modulation index computed as KL/np.log(nbins)

    Example
    -------
    PLV = compute_plv_count(data, 5000.0)

    

=======

# epycom.univariate.compute_pse

    Power spectral entropy

    Parameters
    ----------
    sig: np.array
        time series (float)

    Returns
    -------
    pse - power spectral entropy of analyzed signal, a non-negative value

    Example
    -------
    pac = comute_pse(sig)
    

=======

# epycom.univariate.compute_sample_entropy

       Function to compute sample entropy

       Parameters
       ----------
       sig: np.ndarray
           1D signal
       r: np.float64
           filtering treshold, recommended values: (0.1-0.25)*np.nanstd(sig)
       m: int
           window length of compared run of data, recommended (2-8)

       Returns
       -------
       entropy: numpy.float64 (computed as -np.log(A / B))
           approximate entropy

       Example
       -------
       sample_entropy = approximate_entropy(data, 0.1*np.nanstd(data))
    

=======

# epycom.univariate.compute_shanon_entropy

    Fucntion computes shannon entropy of given signal

    Parameters
    ----------
    sig: np.ndarray
        Signal to analyze

    Returns
    -------
    entro: np.float64
        Computed Shannon entropy of given signal
    

=======

# epycom.univariate.compute_signal_stats

    Function to analyze basic signal statistics

    Parameters:
    ----------
    sig: np.array
        signal to analyze, time series (array, int, float)

    Returns
    -------
    results: tuple
        - power_std: standard deviation of power in band
        - power_mean: mean of power in band
        - power_median: median of power in band
        - power_max: max value of power in band
        - power_min: min value of power in band
        - power_perc25: 25 percentile of power in band
        - power_perc75: 75 percentile of power in band

    Example
    -------
    sig_stats = compute_signal_stats(sig)
    

=======

# epycom.univariate.hjorth_complexity


=======

# epycom.univariate.hjorth_mobility


=======

# epycom.univariate.lyapunov_exponent


=======

# epycom.univariate.mean_vector_length


=======

# epycom.univariate.modulation_index


=======

# epycom.univariate.phase_locking_value


=======

# epycom.univariate.power_spectral_entropy


=======

# epycom.univariate.sample_entropy


=======

# epycom.univariate.shannon_entropy


=======

# epycom.univariate.signal_stats


=======

