function A = inverse_dist(locs)

n = size(locs,1);

A = nan(n,n);
%Abin = nan(n,n);

% Get the distance between adjacent locs
%diff_locs = vecnorm(diff(locs,1,1),2,2);
%median_dist = nanmedian(diff_locs);
%thresh = ceil(median_dist);

for i = 1:n
    for j = 1:i-1
        dist = vecnorm(locs(i,:)-locs(j,:));
        inv_dist = 1/abs(dist).^2;
        
            %{
        if dist < thresh
            Abin(i,j) = 1;
            Abin(j,i) = 1;
        end
            %}
        
        A(i,j) = inv_dist;
        A(j,i) = inv_dist;
        
    end
end

if 0
    figure
    tiledlayout(1,2)
    nexttile
    turn_nans_gray(A)
    
    nexttile
    turn_nans_gray(Abin)
end
  

% Wrap it for storage
A = wrap_or_unwrap_adjacency_fc_toolbox(A);
       

end