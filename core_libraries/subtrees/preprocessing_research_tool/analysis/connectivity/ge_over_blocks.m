function ge = ge_over_blocks(fc)

%{
This takes a 3 dimensional (nchxnchxnblocks) adjacency matrix and outputs
an nblocks x 1 vector of global efficiencies
%}

nblocks = size(fc,3);
ge = nan(nblocks,1);
for b = 1:nblocks
    W = fc(:,:,b);
    ge(b) = efficiency_wei(W);
end

end