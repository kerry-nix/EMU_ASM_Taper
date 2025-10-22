function [mb,lb] = cross_correlation(values,fs,tw,do_tw)

max_lag = 200; % 200 ms
ml = round(max_lag * 1e-3 *fs);

if do_tw

    % divide into time windows
    tw = round(tw*fs);
    times = 1:tw:size(values,1);
    
    % Prep the variables
    mb_all = nan(size(values,2),size(values,2),length(times)-1);
    lb_all = nan(size(values,2),size(values,2),length(times)-1);
    
    for t = 1:length(times)-1
    
        % Calculate cross correlation. 'Normalized' normalizes the sequence so that
        % the autocorrelations at zero lag equal 1
        [r,lags] = xcorr(values(times(t):times(t+1),:),ml,'normalized');
        
        % find the maximum xcorr, and its corresponding lag, for each channel pair
        [M,I] = max(r,[],1);
        nlags = lags(I);
        
        % back out what channels are what
        mb = reshape(M,size(values,2),size(values,2));
        lb = reshape(nlags,size(values,2),size(values,2));
        
        % Make anything that is a nan in mb a nan in lb
        lb(isnan(mb)) = nan;
        
        % make anything at the edge of lags nan in lb
        lb(lb==ml) = nan;
        lb(lb==-ml) = nan;
        
        lb = lb/fs;
    
        mb_all(:,:,t) = mb;
        lb_all(:,:,t) = lb;
    
    end
    
    mb = nanmean(mb_all,3);
    lb = nanmean(lb_all,3);


else
    max_lag = 200; % 200 ms
    ml = round(max_lag * 1e-3 *fs);
    
    % Calculate cross correlation. 'Normalized' normalizes the sequence so that
    % the autocorrelations at zero lag equal 1
    [r,lags] = xcorr(values,ml,'normalized');
    
    % find the maximum xcorr, and its corresponding lag, for each channel pair
    [M,I] = max(r,[],1);
    nlags = lags(I);
    
    % back out what channels are what
    mb = reshape(M,size(values,2),size(values,2));
    lb = reshape(nlags,size(values,2),size(values,2));
    
    % Make anything that is a nan in mb a nan in lb
    lb(isnan(mb)) = nan;
    
    % make anything at the edge of lags nan in lb
    lb(lb==ml) = nan;
    lb(lb==-ml) = nan;
    
    lb = lb/fs;


end

if 0
    figure
    set(gcf,'position',[10 10 1400 400])
    nexttile
    turn_nans_gray(mb)
    nexttile
    turn_nans_gray(lb)
    nexttile
    turn_nans_gray(corr(values))
    
end


end