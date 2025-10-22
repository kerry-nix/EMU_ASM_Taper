function result = movingmean(x, k)
if ndims(x) == 1
    % return np.convolve(np.pad(x, (k - 1, 0), mode='edge'), np.ones(k)/k, mode='valid')
    % return np.convolve(x, np.ones(k)/k, mode=mode)
    result = smoothdata(x, 'movmean', k);
else
    result = zeros(size(x));
    for i = 1:size(x,2)
        result(:, i) = smoothdata(x(:,i), 'movmean', k);
    end
end
