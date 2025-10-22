function [bad,details] = identify_bad_chs(values,fs)

%% Parameters to reject super high variance
tile = 99;
mult = 10;
num_above = 1;
abs_thresh = 5e3;

%% Parameter to reject high 60 Hz
percent_60_hz = 0.7;

%% Parameter to reject electrodes with much higher std than most electrodes
mult_std = 10;

nchs = size(values,2);
chs = 1:nchs;
bad = [];
high_ch = [];
nan_ch = [];
zero_ch = [];
high_var_ch = [];
noisy_ch = [];
all_std = nan(nchs,1);


for i = 1:nchs
    
    bad_ch = 0;
    
    ich = i;
    %ich = 43;
    eeg = values(:,ich);
    bl = nanmedian(eeg);
    
    %% Get channel standard deviation
    all_std(i) = nanstd(eeg);
    
    %% Remove channels with nans in more than half
    if sum(isnan(eeg)) > 0.5*length(eeg)
        bad = [bad;ich];
        nan_ch = [nan_ch;ich];
        continue;
    end
    
    %% Remove channels with zeros in more than half
    if sum(eeg == 0) > 0.5 * length(eeg)
        bad = [bad;ich];
        zero_ch = [zero_ch;ich];
        continue;
    end
    
    %% Remove channels with too many above absolute thresh
    
    if sum(abs(eeg - bl) > abs_thresh) > 10
        bad = [bad;ich];
        bad_ch = 1;
        high_ch = [high_ch;ich];
    end
    
    if 0
    figure
    if bad_ch == 1
        plot(eeg,'r')
    else
        plot(eeg,'k')
    end
    hold on
    plot(xlim,[bl + abs_thresh bl+abs_thresh])
    plot(xlim,[bl - abs_thresh bl-abs_thresh])
    title(sprintf('Percent above threshold: %1.1f%%',100*sum(abs(eeg - bl) > abs_thresh)/length(eeg)));
    pause
    close(gcf)
    end
    
    if bad_ch == 1
        continue;
    end
    
    
    %% Remove channels if there are rare cases of super high variance above baseline (disconnection, moving, popping)
    pct = prctile(eeg,[100-tile tile]);
    thresh = [bl - mult*(bl-pct(1)), bl + mult*(pct(2)-bl)];
    sum_outside = sum(eeg > thresh(2) | eeg < thresh(1));
    if sum_outside >= num_above
        bad_ch = 1;
    end
    
    if 0
        if bad_ch == 1
            plot(eeg,'r');
        else
            plot(eeg,'k');
        end
        title(chLabels{ich})
        hold on
        plot(xlim,[bl bl])
        plot(xlim,[thresh(1) thresh(1)]);
        plot(xlim,[thresh(2) thresh(2)]);
        title(sprintf('Sum outside: %d',sum_outside));
        pause
        hold off
    end
    
    if bad_ch == 1
        bad = [bad;ich];
        high_var_ch = [high_var_ch;ich];
        continue;
    end
    
    %% Remove channels with a lot of 60 Hz noise, suggesting poor impedance
    
    % Calculate fft
    %orig_eeg = orig_values(:,ich);
    %Y = fft(orig_eeg-mean(orig_eeg));
    Y = fft(eeg-nanmean(eeg));
    
    % Get power
    P = abs(Y).^2;
    freqs = linspace(0,fs,length(P)+1);
    freqs = freqs(1:end-1);
    
    % Take first half
    P = P(1:ceil(length(P)/2));
    freqs = freqs(1:ceil(length(freqs)/2));
    
    P_60Hz = sum(P(freqs > 58 & freqs < 62))/sum(P);
    if P_60Hz > percent_60_hz
        bad_ch = 1;
    end
    
    if 0
        figure
        subplot(2,1,1)
        plot(eeg)
        
        subplot(2,1,2)
        spectrogram(eeg,[],[],[],fs);
        title(sprintf('%s Percent 60 Hz power %1.1f',chLabels{ich},P_60Hz*100))
        pause
        close(gcf)
    end
    
    if bad_ch == 1
        bad = [bad;ich];
        noisy_ch = [noisy_ch;ich];
        continue;
    end
    
    
    
end

%% Remove channels for whom the std is much larger than the baseline
median_std = nanmedian(all_std);
higher_std = chs(all_std > mult_std * median_std);
bad_std = higher_std;
bad_std(ismember(bad_std,bad)) = [];
bad = ([bad;bad_std]);
bad_bin = zeros(nchs,1);
bad_bin(bad) = 1;
bad = logical(bad_bin);

details.noisy = noisy_ch;
details.nans = nan_ch;
details.zeros = zero_ch;
details.var = high_var_ch;
details.higher_std = bad_std;
details.high_voltage = high_ch;

end