function out = remove_percent(result,percent)
% result, 24x36x36x2
chans = find(~isnan(squeeze(result(1,:,1,1))));
n_chan = length(find(~isnan(squeeze(result(1,:,1,1)))));
n_remove = round(n_chan*percent);
chans_remove = chans(randperm(n_chan,n_remove));
out = result;
out(:,chans_remove,chans_remove,:) = NaN;
