# epycom.bivariate.Coherence


=======

# epycom.bivariate.LinearCorrelation


=======

# epycom.bivariate.PhaseConsistency


=======

# epycom.bivariate.PhaseLagIndex


=======

# epycom.bivariate.PhaseSynchrony


=======

# epycom.bivariate.RelativeEntropy


=======

# epycom.bivariate.SpectraMultiplication


=======

# epycom.bivariate.coherence


=======

# epycom.bivariate.compute_coherence

    Magnitude squared coherence between two time series (raw,
    not filtered signals)

    When win and win_step is not 0, calculates evolution of coherence

    When win>len(sig) or win<=0, calculates only one coherence value

    When lag and lag_step is not 0, shifts the sig[1] from negative
    to positive lag and takes the max coherence (best fit)

    Parameters
    ----------
    sig: np.array
        2D numpy array of shape (signals, samples), time series (int, float)
    fs: int, float
        sampling frequency in Hz
    fband: list
        frequency range in Hz (float)
    lag: int
        negative and positive shift of time series in samples
    lag_step: int
        step of shift in samples
    fft_win: int
        length of fft window in sec

    Returns
    -------
    max_coh: list
        maximum coherence in shift
    tau: float
        shift of maximum coherence in samples,
        value in range <-lag,+lag> (float)
        tau<0: sig[1] -> sig[0]
        tau>0: sig[0] -> sig[1]

    Example
    -------
    max_coh,tau = compute_coherence(sig, fs=5000, fband=[1.0,4.0], lag=0,
                                    lag_step=0, win=0, win_step=0, fft_win=1)
    

=======

# epycom.bivariate.compute_lincorr

    Linear correlation (Pearson's coefficient) between two time series

    When lag and lag_step is not 0, shifts the sig[1] from negative
    to positive lag and takes the max correlation (best fit)

    Parameters
    ----------
    sig: np.array
        2D numpy array of shape (signals, samples), time series (int, float)
    lag: int
        negative and positive shift of time series in samples
    lag_step: int
        step of shift

    Returns
    -------
    lincorr: list
        maximum linear correlation in shift
    tau: float
        shift of maximum correlation in samples,
        value in range <-lag,+lag> (float)
        tau<0: sig[1] -> sig[0]
        tau>0: sig[0] -> sig[1]

    Example
    -------
    lincorr,tau = compute_lincorr(sig, 200, 20)
    

=======

# epycom.bivariate.compute_phase_const


    **under development**

    calculation of phase consistency between two signals
    irrespective of the amplitude
    pre-filtering of the signals is necessary
    use appropriate lag and step (it calculates phase_const between single
    lag steps in whole length of given time signals)

    Parameters
    ----------
    sig: np.array
        2D numpy array of shape (signals, samples), time series (float)
    lag: int
        negative and positive shift of time series in samples
    lag_step: int
        step of shift in samples

    Returns
    -------
    phase_const: float
        ranges between 0 and 1
        (1 for the phase lock which does not shift during the time period)

    Example
    -------
    phsc = compute_phase_const(sig, 500, 50)
    

=======

# epycom.bivariate.compute_phase_sync

    Calculation of phase synchronization using Hilbert transformation
    sensitive to phases, irrespective of the amplitude
    and phase shift, pre-filtering of the signals is necessary

    Parameters
    ----------
    sig: np.array
        2D numpy array of shape (signals, samples), time series ( float)

    Returns
    -------
    phase_sync: float
        ranges between 0 and 1 (1 for the perfect synchronization)

    Example
    -------
    phs = compute_phase_sync(sig)

    References
    ----------
    Quiroga et al. 2008
    

=======

# epycom.bivariate.compute_pli

    Phase-lag index.

    - filter signal before pli calculation (if one is filtered and the other
      is not (or in different f-band), it can return fake high pli, which is
      caused by substraction (difference) of inst. phases at different scales)
    - use appropriate win and lag (max lag ~= fs/2*fmax, else it runs over one
      period and finds pli=1)

    Parameters
    ----------
    sig: np.array
        2D numpy array of shape (signals, samples), time series (float)
    lag: int
        negative and positive shift of time series in samples
    lag_step: int
        step of shift in samples

    Returns
    -------
    pli: float
        ranges between 0 and 1 (1 for the best phase match between signals)
    tau: int
        phase lag for max pli value (in samples, 0 means no lag)

    Example
    -------
    pli, tau = compute_pli(sig, lag=500, lag_step=50)

    References
    ----------
    [1] C. J. Stam and J. C. Reijneveld, “Graph theoretical analysis of
    complex networks in the brain,” Nonlinear Biomed. Phys., vol. 1, no. 1,
    p. 3, 2007.
    

=======

# epycom.bivariate.compute_relative_entropy

    Calculation of Kullback-Leibler divergence:
    relative entropy of sig[0] with respect to sig[1]
    and relative entropy of sig[1] with respect to sig[0]

    Parameters
    ----------
    sig: np.array
        2D numpy array of shape (signals, samples), time series (int, float)

    Returns
    -------
    ren: float
        max value of relative entropy between sig[0] and sig[1]

    Example:
    -------
    ren = compute_relative_entropy(sig)
    

=======

# epycom.bivariate.compute_spect_multp

    Multiply spectra of two time series and transforms it back to time domain,
    where the mean and std is calculated

    Parameters
    ----------
    sig: np.array
        2D numpy array of shape (signals, samples), time series (int, float)

    Returns
    -------
    sig_sm_mean: float
        aritmetic mean value of multiplied signals
    sig_sm_std: float
        standard deviation of multiplied signals

    Example
    -------
    mspect = spect_multp(sig)
    

=======

# epycom.bivariate.linear_correlation


=======

# epycom.bivariate.phase_consistency


=======

# epycom.bivariate.phase_lag_index


=======

# epycom.bivariate.phase_synchrony


=======

# epycom.bivariate.relative_entropy


=======

# epycom.bivariate.spectra_multiplication


=======

