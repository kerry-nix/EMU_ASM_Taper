function all_plv = plv_calc(values,fs,tw,do_tw)

%% Get filtered signal
out = filter_canonical_freqs(values,fs); % broadband at end
nfreqs = size(out,3);
nchs = size(values,2);


if do_tw
    tw = round(tw*fs);
    times = 1:tw:size(values,1);
    all_plv = nan(nchs,nchs,nfreqs,length(times)-1);
    for t = 1:length(times)-1
        for f = 1:nfreqs
            filteredData = out(times(t):times(t+1),:,f);
            % Get phase of each signal
            phase = nan(size(filteredData));
            for ich = 1:nchs
                phase(:,ich)= angle(hilbert(filteredData(:,ich)));
            end
            % Get PLV
            plv = nan(nchs,nchs);
            for ich = 1:nchs
                for jch = ich+1:nchs
                    e = exp(1i*(phase(:,ich) - phase(:,jch)));
                    plv(ich,jch) = abs(sum(e,1))/size(phase,1);
                    plv(jch,ich) = abs(sum(e,1))/size(phase,1);
                end
            end
            all_plv(:,:,f,t) = plv;

        end
    end
    all_plv = nanmean(all_plv,4);

else
    %% initialize output vector
    all_plv = nan(nchs,nchs,nfreqs);
    
    
    
    % Do plv for each freq
    for f = 1:nfreqs
        
        filteredData = out(:,:,f);
    
        % Get phase of each signal
        phase = nan(size(filteredData));
        for ich = 1:nchs
            phase(:,ich)= angle(hilbert(filteredData(:,ich)));
        end
    
        % Get PLV
        plv = nan(nchs,nchs);
        for ich = 1:nchs
            for jch = ich+1:nchs
                e = exp(1i*(phase(:,ich) - phase(:,jch)));
                plv(ich,jch) = abs(sum(e,1))/size(phase,1);
                plv(jch,ich) = abs(sum(e,1))/size(phase,1);
            end
        end
        %temp_plv(:,:,f) = plv;
        all_plv(:,:,f) = plv;
    
        if 0
            figure
            nexttile
            plot(filteredData(:,2))
            hold on
            plot(filteredData(:,3))
            nexttile
            plot((hilbert(filteredData(:,2))))
            hold on
            plot(hilbert(filteredData(:,3)))
    
            nexttile
            plot(phase(:,2))
            hold on
            plot(phase(:,3))
        end
    end
end



if 0
    figure; tiledlayout(1,6)
    for i = 1:6
        nexttile
        turn_nans_gray(all_plv(:,:,i))
        colorbar
    end
end

end