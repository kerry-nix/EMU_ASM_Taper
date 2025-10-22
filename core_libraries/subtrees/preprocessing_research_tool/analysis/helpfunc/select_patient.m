function select_patient(dataset)

load(which(strcat(dataset,'_metaDataFull.mat')))
% for every file
for i = 1:length(metaData)

    % clean channels
    [chans,~,~] = clean_labels(metaData(i).channels);

    % join to a string
    chanstr = strjoin(chans);
    judge = []; 
    judgemat = zeros(3,12); % electrodes present for A, B, C
    tofind = {metaData(i).tofindA, metaData(i).tofindB, metaData(i).tofindC};
    for ind = 1:3
        for j = 1:12
            if ismember(strcat('L',tofind{ind},num2str(j)),chans) && ismember(strcat('R',tofind{ind},num2str(j)),chans) 
                judgemat(ind,j) = 1;
            end
        end
    end

    if any(any(judgemat))
        metaData(i).LR = 1;  
    else
        metaData(i).LR = 0;
    end
    metaData(i).elecInd = judgemat;
end

% remove files without LR paired chans
subMeta = metaData;
count = 1;
while count <= length(subMeta)
    if subMeta(count).LR == 0
        subMeta(count) = [];
    else
        count = count + 1;
    end
end

% count # of available electrodes and which is full in ABC 
for i = 1:length(subMeta)
    subMeta(i).count = length(find(subMeta(i).elecInd));
    ABC = 'ABC';
    subMeta(i).ABC = [];
    for j = 1:3
        if length(find(subMeta(i).elecInd(j,:))) == 12
            subMeta(i).ABC = [subMeta(i).ABC, ABC(j)];
        end
    end
end
patientList = subMeta;
save(strcat('meta/',dataset,'_patientList.mat'),'patientList')
save(strcat('meta/',dataset,'_subMeta.mat'),'subMeta')

% now all other info were manually added
% load soz info
% soz = readtable('soz.csv');
% soz = table2cell(soz);
% soz(:,2) = cellfun(@(x) lower(x),soz(:,2),'UniformOutput',false);
% for i = 1:length(subMeta)
%     tmp = strsplit(subMeta(i).filename,'_'); 
%     name = tmp{1};
%     [~,ind] = ismember(name,soz(:,1));
%     if ind > 0
%         subMeta(i).soz = soz{ind,2};
%     end
% end

% select files with full 36 chans
% subSubMeta = subMeta;
% count = 1;
% while count <= length(subSubMeta)
%     if ~strcmp(subSubMeta(count).ABC,'ABC')
%         subSubMeta(count) = [];
%     else
%         count = count + 1;
%     end
% end
