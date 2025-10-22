function out = filter_canonical_freqs(values,fs)
% This gives broadband back at the end
%% Parameters
order = 4;
freqs = get_frequencies; 
freqs = freqs(1:end-1,:); % remove the broadband
nchs = size(values,2);

for ich = 1:nchs
    curr_values = values(:,ich);
    curr_values(isnan(curr_values)) = nanmean(curr_values);
    values(:,ich) = curr_values;
end

nfreqs = size(freqs,1);
out = nan(size(values,1),size(values,2),nfreqs+1);
for f = 1:nfreqs
    out(:,:,f) = bandpass_any(values,fs,[freqs(f,1),min(fs/2-1,freqs(f,2))],order);
end

% Put the raw signal back at the END
out(:,:,end) = values;


L = size(values,1);
fr = fs*(0:(L/2))/L;

if 0
    figure
    ch = 7;
    set(gcf,'position',[10 10 1400 400])
    t = tiledlayout(2,nfreqs+1,'TileSpacing','tight','padding','tight');
    nexttile 
    plot(values(:,ch))
    for f = 1:nfreqs
        nexttile
        plot(out(:,ch,f))
        title(sprintf('%d-%d Hz',freqs(f,1),freqs(f,2)))
       
    end
    nexttile
    Y = fft(values(:,ch));
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    plot(fr(fr<70),P1(fr<70)) 
    hold on
    for f = 1:nfreqs
        nexttile
        Y = fft(out(:,ch,f));
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        plot(fr(fr<70),P1(fr<70)) 
        hold on
    end
end

end