function values = notch_filter(values,fs)

f = designfilt('bandstopiir','FilterOrder',4, ...
   'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
   'DesignMethod','butter','SampleRate',fs);
%fvtool(f)
for i = 1:size(values,2) 
    
    eeg = values(:,i);
    
    if sum(~isnan(eeg)) == 0
        continue
    end
    
    eeg(isnan(eeg)) = nanmean(eeg);
    eeg = filtfilt(f,eeg);   
    values(:,i) = eeg;
end


end

