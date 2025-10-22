function out = calcConn(data,fs,method,ind,tw,do_tw)

out = NaN(36,36,1,2);
toFill = find(ind');
leftData = data(:,1:size(data,2)/2);
rightData = data(:,size(data,2)/2+1:end);
switch method
    case 'pearson'
        out(toFill,toFill,:,1) = new_pearson_calc(leftData,fs,tw,do_tw);
        out(toFill,toFill,:,2) = new_pearson_calc(rightData,fs,tw,do_tw);
    case 'squaredPearson'
        out(toFill,toFill,:,1) = new_pearson_calc(leftData,fs,tw,do_tw);
        out(toFill,toFill,:,2) = new_pearson_calc(rightData,fs,tw,do_tw);
        out = out.^2;
    case 'crossCorr'
        out(toFill,toFill,:,1) = cross_correlation(leftData,fs,tw,do_tw);
        out(toFill,toFill,:,2) = cross_correlation(rightData,fs,tw,do_tw);
    case 'coh'
        out = NaN(36,36,7,2);
        out(toFill,toFill,:,1) = faster_coherence_calc(leftData,fs,tw,do_tw);
        out(toFill,toFill,:,2) = faster_coherence_calc(rightData,fs,tw,do_tw);
    case 'plv'
        out = NaN(36,36,7,2);
        out(toFill,toFill,:,1) = plv_calc(leftData,fs,tw,do_tw);
        out(toFill,toFill,:,2) = plv_calc(rightData,fs,tw,do_tw);
    case 'relaEntropy'
        out = NaN(36,36,7,2);
        out(toFill,toFill,:,1) = relative_entropy(leftData,fs,tw,do_tw);
        out(toFill,toFill,:,2) = relative_entropy(rightData,fs,tw,do_tw);
    case 'granger'
        out(toFill,toFill,:,1) = granger(leftData);
        out(toFill,toFill,:,2) = granger(rightData);
end
