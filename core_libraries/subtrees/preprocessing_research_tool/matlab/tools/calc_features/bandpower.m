function bp = bandpower(data, fs, band,varargin)
% Adapted from https://raphaelvallat.com/bandpower.html
% Compute the average power of the signal x in a specific frequency band.
%
%     Parameters
%     ----------
%     data : 1d-array or 2d-array
%         Input signal in the time-domain. (time by channels)
%     fs : float
%         Sampling frequency of the data.
%     band : list
%         Lower and upper frequencies of the band of interest.
%     win : float
%         Length of each window in seconds.
%         If None, win = (1 / min(band)) * 2
%     relative : boolean
%         If True, return the relative power (= divided by the total power of the signal).
%         If False (default), return the absolute power.
%
%     Return
%     ------
%     bp : float
%         Absolute or relative band power.
low = band(1); high = band(2);
win = (1 / low) * 2;
relative = false;
para = varargin;
assert(mod(length(para),2)==0,'CNTtools:invalidInput','Optional inputs should be paired.');
while length(para) >= 2
    name = para{1};
    val = para{2};
    para = para(3:end);
    switch name
        case 'win'
            win = val;
        case 'relative'
            relative = val;
    end
end


% Define window length
nperseg = win * fs;

% Compute the modified periodogram (Welch)
[pxx, f] = pwelch(data, nperseg, [], [], fs);
% Frequency resolution
freq_res = f(2) - f(1);

% Find closest indices of band in frequency vector
idx_band = find(f >= low & f <= high);

% Integral approximation of the spectrum using Simpson's rule.
if ismatrix(pxx)
    bp = simpson(freq_res,pxx(idx_band,:));
else
    bp = simpson(freq_res,pxx(idx_band));
end
if relative
    bp = bp/simpson(freq_res,pxx);
end

