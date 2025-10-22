function re = relative_entropy(values,fs,tw,do_tw)

%% Get filtered signal
out = filter_canonical_freqs(values,fs); % broadband is at the end
nfreqs = size(out,3);
nchs = size(values,2);


if do_tw
    % divide into time windows
    tw = round(tw*fs);
    times = 1:tw:size(values,1);
    re = nan(nchs,nchs,nfreqs,length(times)-1);
    
    % loop over times and frequences
    for t = 1:length(times)-1
        for f = 1:nfreqs
        
            filteredData = out(times(t):times(t+1),:,f);
        
            for ich = 1:nchs
                for jch = ich+1:nchs
                    %h1 = histcounts(filteredData(:,ich),10);
                    %h2 = histcounts(filteredData(:,jch),10);

                    h1 = (steve_histcounts(filteredData(:,ich),10))'; % faster
                    h2 = (steve_histcounts(filteredData(:,jch),10))';
                    smooth = 1e-10;
                    h1 = h1 + smooth; h2 = h2 + smooth;
                    h1 = h1/sum(h1); h2 = h2/sum(h2); % normalize?
                    S1 = sum(h1.*log(h1./h2));
                    S2 = sum(h2.*log(h2./h1));
                    re(ich,jch,f,t) = max([S1,S2]);
                    re(jch,ich,f,t) = re(ich,jch,f,t);
                end
            end
        end
    end
    
    re = nanmean(re,4);

else
    re = nan(nchs,nchs,nfreqs);
    
    
    % loop over frequences
    for f = 1:nfreqs
    
        filteredData = out(:,:,f);
    
        for ich = 1:nchs
            for jch = ich+1:nchs
                h1 = (steve_histcounts(filteredData(:,ich),10))';
                h2 = (steve_histcounts(filteredData(:,jch),10))';

                %h1 = histcounts(filteredData(:,ich),10);
                %h2 = histcounts(filteredData(:,jch),10);
                h1 = h1/sum(h1); h2 = h2/sum(h2); % normalize?
                %S1 = sum(h1*log(h1/h2));
                %S2 = sum(h2*log(h2/h1));
                S1 = sum(h1.*log(h1./h2));
                S2 = sum(h2.*log(h2./h1));
                re(ich,jch,f) = max([S1,S2]);
                re(jch,ich,f) = re(ich,jch,f);
            end
        end
    end
end

if 0
    figure
    for i = 1:nfreqs
    nexttile
    turn_nans_gray(re(:,:,i))
    end

end
%}


end