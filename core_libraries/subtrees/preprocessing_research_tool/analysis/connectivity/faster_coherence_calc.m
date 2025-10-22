function all_coherence = faster_coherence_calc(values,fs,tw,do_tw)

%% Parameters
window = fs * 1;
NFFT = window;
freqs = get_frequencies; 
nfreqs = size(freqs,1);

nchs = size(values,2);

%% initialize output vector
all_coherence = nan(nchs,nchs,nfreqs);

for ich = 1:nchs
    curr_values = values(:,ich);
    curr_values(isnan(curr_values)) = nanmean(curr_values);
    values(:,ich) = curr_values;
end


%% Remove nan rows, keeping track
nan_rows = any(isnan(values),1); % find channels with nans for any time points
values_no_nans = values(:,~nan_rows);
nchs_no_nans = size(values_no_nans,2);


if do_tw
    tw = round(tw*fs);
    times = 1:tw:size(values,1);
    temp_coherence = nan(nchs_no_nans,nchs_no_nans,nfreqs,length(times)-1);
    for t = 1:length(times)-1
        for ich = 1:nchs_no_nans
            [cxy,f] = mscohere(values_no_nans(times(t):times(t+1),ich),values_no_nans(times(t):times(t+1),:),hamming(window),[],NFFT,fs);
            for i_f = 1:nfreqs
                temp_coherence(:,ich,i_f,t) = ...
                    nanmean(cxy(f >= freqs(i_f,1) & f <= freqs(i_f,2),:),1);
            end
        end

    end
    temp_coherence = nanmean(temp_coherence,4);
else
    temp_coherence = nan(nchs_no_nans,nchs_no_nans,nfreqs);
    for ich = 1:nchs_no_nans

        % Do MS cohere on full thing
        [cxy,f] = mscohere(values_no_nans(:,ich),values_no_nans,hamming(window),[],NFFT,fs);
    
        % Average coherence in frequency bins of interest
        for i_f = 1:nfreqs
            temp_coherence(:,ich,i_f) = ...
                nanmean(cxy(f >= freqs(i_f,1) & f <= freqs(i_f,2),:),1);
        end

    end
    
end

%% Put the non-nans back
all_coherence(~nan_rows,~nan_rows,:) = temp_coherence;
all_coherence(logical(repmat(eye(nchs,nchs),1,1,nfreqs))) = nan;


end