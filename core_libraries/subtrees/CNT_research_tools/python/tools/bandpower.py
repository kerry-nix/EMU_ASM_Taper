
from scipy.signal import welch
from scipy.integrate import simpson
import numpy as np


def bandpower(data, sf, band, window_sec=None, relative=False):
    """Adapted from https://raphaelvallat.com/bandpower.html
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
    """
    band = np.asarray(band)
    low, high = band

    # Define window length
    if window_sec is not None:
        nperseg = window_sec * sf
    else:
        nperseg = (2 / low) * sf

    # Compute the modified periodogram (Welch)
    freqs, psd = welch(data.T, sf, nperseg=nperseg)

    # Frequency resolution
    freq_res = freqs[1] - freqs[0]

    # Find closest indices of band in frequency vector
    idx_band = np.logical_and(freqs >= low, freqs <= high)

    # Integral approximation of the spectrum using Simpson's rule.
    if psd.ndim == 2:
        bp = simpson(psd[:, idx_band], dx=freq_res)
    elif psd.ndim == 1:
        bp = simpson(psd[idx_band], dx=freq_res)

    if relative:
        bp /= simpson(psd, dx=freq_res)
    return bp