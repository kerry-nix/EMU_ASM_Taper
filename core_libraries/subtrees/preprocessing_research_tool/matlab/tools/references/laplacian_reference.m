function [out_values,close_chs,laplacian_labels] = laplacian_reference(values,include,locs,radius,labels)

assert(~isempty(locs))

out_values = nan(size(values));
nchs = size(values,2);

D = pdist2(locs,locs); 
close = D < radius;
close(logical(eye(size(close)))) = 0;

close_chs = cell(nchs,1);

for i = 1:nchs
    out_values(:,i) = values(:,i) - nanmean(values(:,close(i,:)&include'),2);
    close_chs{i} = find(close(i,:)&include');
end

%laplacian_labels = strcat(labels,'-laplacian');
laplacian_labels = cell(length(labels),1);
for i = 1:length(laplacian_labels)
    ch_string = '';
    for j = 1:length(close_chs{i})
        ch_string = [ch_string,labels{close_chs{i}(j)}];
        if j < length(close_chs{i})
            ch_string = [ch_string,', '];
        end
    end
    laplacian_labels{i} = [labels{i},'-',ch_string];
end


end