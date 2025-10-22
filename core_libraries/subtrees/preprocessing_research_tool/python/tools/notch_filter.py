from scipy.signal import iirnotch, filtfilt, butter



def notch_filter(data,fs):

    b, a = iirnotch(60, 30, fs=fs)
    # b, a = butter(4, [59,61], 'bandstop', fs) 

    y = filtfilt(b, a, data, axis = 0)

    return y

