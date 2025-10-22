function target_thing = reconcile_ch_names(target_chs,new_chs,new_thing)

% clean both
target_chs = clean_labels(target_chs);
new_chs = clean_labels(new_chs);

[Lia,locb] = ismember(target_chs,new_chs);
target_thing = nan(size(target_chs,1),size(new_thing,2));
target_thing(Lia,:) = new_thing(locb(Lia),:);
expanded_chs = cell(length(target_chs),1);
expanded_chs(Lia) = new_chs(locb(Lia));

if 0
   T= table(target_chs,expanded_chs,target_thing)
end

end