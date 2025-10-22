import numpy as np
from beartype import beartype
from beartype.typing import Iterable

@beartype
def common_average_ref(data: np.ndarray, labels: Iterable[str]):
    """Do common average reference for the input data array

    Args:
        data (np.ndarray): _description_
        fs (Number): _description_

    Returns:
        list,dict: _description_
    """
    out_data = data - np.nanmean(data,1)[:,np.newaxis]
    car_labels = [label+'-CAR' for label in labels]

    return out_data,car_labels