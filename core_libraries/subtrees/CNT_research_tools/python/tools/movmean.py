import numpy as np
from scipy.ndimage.filters import uniform_filter1d

def movmean(x, k):
    if x.ndim == 1:
        # return np.convolve(np.pad(x, (k - 1, 0), mode='edge'), np.ones(k)/k, mode='valid')
        # return np.convolve(x, np.ones(k)/k, mode=mode)
        return uniform_filter1d(x, size=k)
    else:
        avgd_x = np.zeros(x.shape)
        for i, row in enumerate(x):
            # avgd_x[i, :] = np.convolve(row, np.ones(k)/k, mode=mode)
            avgd_x[i, :] = uniform_filter1d(row, size=k)
        return avgd_x
if __name__ == "__main__":
    print(movmean(np.random.rand(100), k=10))