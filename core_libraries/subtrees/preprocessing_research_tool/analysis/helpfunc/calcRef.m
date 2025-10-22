function [values,labels,newInd,keep] = calcRef(values,labels,reference,oldInd)

%% Common average reference (include only intra-cranial)
switch reference
    case 'car'
        [values,labels] = common_average_reference(values,true(length(labels),1),labels);
        labels = cellfun(@(x) x(1:regexp(x,'-')-1),labels,'UniformOutput',false);
        newInd = oldInd;
        keep = 1;
    case 'bipolar'
        [values,~,newlabels,~] = bipolar_montage(values,labels);
        ind = cellfun(@(x) strcmp(x,'-'),newlabels);
        % newlabels = cellfun(@(x) x(1:regexp(x,'-')-1),labels,'UniformOutput',false);
        [values,labels,newInd,keep] = updateLabel(values,labels,ind,oldInd);
        % labels = labels(cellfun(@(x) ismember(x(1:regexp(x,'-')-1),newlabels),labels));
end


