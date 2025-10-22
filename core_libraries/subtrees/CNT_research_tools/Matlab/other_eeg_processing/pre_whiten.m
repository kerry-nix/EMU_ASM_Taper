function data = pre_whiten(data)


for j = 1:size(data,2)
    vals = data(:,j);    
    
    if sum(~isnan(vals)) == 0
        continue
    end
    
    %{
    mdl = regARIMA('ARLags',1);
    mdl = estimate(mdl,vals,'Display','Off');
    E = infer(mdl,vals);
    %}

    %
    mdl = fitlm(vals(1:end-1),vals(2:end));
    E = mdl.predict(vals(1:end-1)) - vals(2:end);
    %}

    if length(E) < length(vals)
        E = [E;nan(length(vals)-length(E),1)];
    end
    data(:,j) = E;
end



end