function [values,newlabels,newInd,keep] = updateLabel(values,labels,bad,oldInd)
% this function updates the 3x12 ind mat, if keep the patient, newlabels,
% values
% load('elecmat.mat');
% 1. find bad channels, set the ind in elecInd to 0
num_chan = length(labels);
keep = 1;
% update labels
badLeft = bad(1:num_chan/2);
badRight = bad(num_chan/2+1:end);
badFull = badLeft | badRight;
badFull = [badFull;badFull];
newlabels = labels(~badFull);
% update elecInd
bad_chan = labels(bad);
[bad_chan, elec, num] = clean_labels(bad_chan);
newInd = oldInd;
for i = 1:length(bad_chan)
    elecName = elec{i};
    elecName = elecName(2:end);
    if strcmp(elecName,'DA')|strcmp(elecName,'AD')|strcmp(elecName,'A')
        newInd(1,num(i)) = 0;
    elseif strcmp(elecName,'DH')|strcmp(elecName,'HD')|strcmp(elecName,'B')|strcmp(elecName,'AH')
        newInd(2,num(i)) = 0;
    elseif strcmp(elecName,'C')|strcmp(elecName,'PH')
        newInd(3,num(i)) = 0;
    end
end
% update values
values = values(:,ismember(labels,newlabels));
if isempty(newlabels)
    keep = 0;
end