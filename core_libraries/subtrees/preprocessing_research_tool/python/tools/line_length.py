import numpy as np


def line_length(signal: np.ndarray):
    """_summary_

    Args:
        signal (np.ndarray): _description_

    Returns:
        _type_: _description_
    """
    return np.mean(np.abs(np.diff(signal, axis=0)), axis=0)
