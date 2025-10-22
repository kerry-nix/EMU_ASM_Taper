function out = bandpass_any(values,fs,range,order)

d = designfilt('bandpassiir','FilterOrder',order, ...
    'HalfPowerFrequency1',range(1),'HalfPowerFrequency2',range(2), ...
    'SampleRate',fs);
out = nan(size(values));

for i = 1:size(values,2)
    eeg = values(:,i);
    
    if sum(~isnan(eeg)) == 0
        continue
    end
    
    eeg(isnan(eeg)) = nanmean(eeg); % make nans equal to the mean.
    eeg = filtfilt(d,eeg);   
    out(:,i) = eeg;
end

end