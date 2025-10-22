function c = steve_histcounts(x, n)
    L = min(x);
    U = max(x);
    W = (U-L)/n;
    i = min(n, 1+floor((x-L)/W));
    c = accumarray(i, 1, [n 1]);
end