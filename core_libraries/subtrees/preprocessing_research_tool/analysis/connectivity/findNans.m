function [judgeNan,percentNan] = findNans(data)
judgeNan = 0;
percentNan = 0;
nanMat = isnan(data);
if any(sum(nanMat,1) > 0.2*size(data,1))
    judgeNan = 1;
    percentNan = sum(sum(nanMat,1) > 0.2*size(data,1))/size(data,2);
end