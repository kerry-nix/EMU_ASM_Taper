function ns = ns_over_blocks(fc)

%{
This takes a 3 dimensional (nchxnchxnblocks) adjacency matrix and outputs
an nblocks x 1 vector of global efficiencies
%}

nblocks = size(fc,3);
nchs = size(fc,1);
ns = nan(nchs,nblocks);
for b = 1:nblocks
    W = fc(:,:,b);
    W(logical(eye(size(W)))) = nan;
    nan_rows = isnan(mean(W,1,'omitnan'));
    ns_int = sum(W,1,'omitnan');
    ns_int(nan_rows) = nan;
    ns(:,b) = ns_int;
    
end

end