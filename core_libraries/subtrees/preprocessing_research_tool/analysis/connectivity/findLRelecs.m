function [judgemat,keep] = findLRelecs(chans)
chanstr = strjoin(chans);
judgemat = zeros(3,12);
elecInd = char(65:67);
for ind = 1:3
    for j = 1:12
        if ismember(strcat('L',elecInd(ind),num2str(j)),chans) && ismember(strcat('R',elecInd(ind),num2str(j)),chans) 
            judgemat(ind,j) = 1;
        end
    end
end
if any(any(judgemat))
    keep = 1;
else
    keep = 0;
end
