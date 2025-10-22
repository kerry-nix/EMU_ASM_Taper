function all_coherence = coherence_calc(values,fs)

%% Parameters
hws = 2;
overlaps = 1;
freqs = get_frequencies;
nfreqs = size(freqs,1);

nchs = size(values,2);

% convert hamming window and overlap from seconds to segments
hw = round(hws*fs);
overlap = round(overlaps*fs);

%% initialize output vector
all_coherence = nan(nchs,nchs,nfreqs);

for ich = 1:nchs
    curr_values = values(:,ich);
    curr_values(isnan(curr_values)) = nanmean(curr_values);
    values(:,ich) = curr_values;
end

for ich = 1:nchs
    
    if any(isnan(values(:,ich))), continue; end
    
    for jch = 1:ich-1
        
        if any(isnan(values(:,jch))), continue; end
        
        %% Calculate coherence
        [cxy,f] = mscohere(values(:,ich),values(:,jch),hw,overlap,[],fs);
        
        %% Average coherence in frequency bins of interest
        for i_f = 1:nfreqs
            all_coherence(ich,jch,i_f) = ...
                nanmean(cxy(f >= freqs(i_f,1) & f <= freqs(i_f,2)));
            
            all_coherence(jch,ich,i_f) = ...
                nanmean(cxy(f >= freqs(i_f,1) & f <= freqs(i_f,2)));
        end
        
    end
end


end