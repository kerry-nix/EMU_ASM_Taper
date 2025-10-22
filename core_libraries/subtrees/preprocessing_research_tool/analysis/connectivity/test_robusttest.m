% try redo reference for 1 patient
ind = find(triu(ones(36,36), 1));
perc = [0.2 0.4 0.6 0.8];
n_trial = 100;
result = NaN(length(refMethods),sum(numFeats),36,36,2);
permResult = NaN(48,length(perc),n_trial,72);%for each patient X 48 methods X 3 percent X 500 trials  X n contacts
globalPermResult = NaN(48,length(perc),n_trial);%for each method X 3 percent X 500 trials X 62 patients;
filename = strcat('data/',patientList(46).patient,'.mat');
for i = 1:4
    disp(i)
    for n = 1:n_trial
        load(filename)
        n_chan = length(labels)/2;
        n_remove = round(n_chan*perc(i)/2);
        remove = randperm(n_chan,n_remove);
        labels([remove,remove+n_chan]) = [];
        data(:,[remove,remove+n_chan]) = [];
        for j = 1:length(refMethods)
            [values,newlabels,elecIndNew,keepPatient] = calcRef(data,labels,refMethods{j});
            if keepPatient == 0
                continue
            end
            for k = 1:length(connMeasures)
                out = calcConn(values,fs,connMeasures{k},elecIndNew,2,1);
                result(j,sum(numFeats(1:k-1))+1:sum(numFeats(1:k-1))+numFeats(k),:,:,:) = permute(out,[3,1,2,4]);
            end
        end
        tmp = reshape(result,[length(refMethods)*sum(numFeats),36,36,2]);
        nodeStr = [];
        for j = 1:36
            nodeStr(:,j,1) = mean(tmp(:,j,[1:(j-1),(j+1):end],1),3,'omitnan');
            nodeStr(:,j,2) = mean(tmp(:,j,[1:(j-1),(j+1):end],2),3,'omitnan');
        end
        nodeStr = reshape(nodeStr,[length(refMethods)*sum(numFeats),36*2]);
        globalStr = mean(nodeStr,2,'omitnan');
        permResult(:,i,n,:) = nodeStr;
        globalPermResult(:,i,n) = globalStr;
    end
end
permReliability = NaN(48,length(perc));% for each patient X 48 method X 3 percent has one reliability score;
% This is the variance across random permutations, averaged over all
% electrodes
var_error = nanmean(squeeze(nanstd(permResult,0,3).^2),3);
% This is the variance across electrodes for a given permutation of the
% graph at that removal percentage, then averaged over all permutations for
% that removal percentage
var_true = nanmean((nanstd(permResult,0,4).^2),3);
permReliability = var_true./(var_true + var_error);
save('data/test_perm_reref.mat','permResult','globalPermResult','var_error','var_true','permReliability')
% original method
load(which('results.mat'))
permResult = NaN(48,length(perc),n_trial,72);%for each patient X 48 methods X 3 percent X 500 trials  X n contacts
for k = 1:4
    for n = 1:n_trial
        tmp = results(30,:,:,:,:,:);
        tmp = reshape(tmp,[length(refMethods),sum(numFeats),36,36,2]);
        for j = 1:2
            tmp(j,:,:,:,:) = remove_percent(squeeze(tmp(j,:,:,:,:)),perc(k));
        end
        tmp = permute(tmp,[2,1,3,4,5]);
        tmp = reshape(tmp,[length(refMethods)*sum(numFeats),36,36,2]);
        nodeStr = [];
        for j = 1:36
            nodeStr(:,j,1) = mean(tmp(:,j,[1:(j-1),(j+1):end],1),3,'omitnan');
            nodeStr(:,j,2) = mean(tmp(:,j,[1:(j-1),(j+1):end],2),3,'omitnan');
        end
        nodeStr = reshape(nodeStr,[length(refMethods)*sum(numFeats),36*2]);
        globalStr = mean(nodeStr,2,'omitnan');
        permResult(:,k,n,:) = nodeStr;
    end
end
% save('data/perm.mat','permResult','globalPermResult')
permReliability_org = NaN(48,length(perc));% for each patient X 48 method X 3 percent has one reliability score;
% This is the variance across random permutations, averaged over all
% electrodes
var_error  = nanmean(squeeze(nanstd(permResult,0,3).^2),3);
% This is the variance across electrodes for a given permutation of the
% graph at that removal percentage, then averaged over all permutations for
% that removal percentage
var_true  = nanmean((nanstd(permResult,0,4).^2),3);
permReliability_org  = var_true./(var_true + var_error);
save('data/test_perm_reref.mat','permResult','var_error','var_true','permReliability','permReliability_org')
%% add error bar, beautiful/consistent color, etc.
subplot(1,2,1)
plot(permReliability)
legend('20','40','60','80','FontSize',14)
subplot(1,2,2)
plot(permReliability_org)
legend('20','40','60','80')
%% SPLIT BY REF
ref = reshape(permReliability,[63,24,2,4]);
ref = squeeze(nanmean(ref,2));
blue = hex2rgb('#51b8bd');red = hex2rgb('#de7862');
col = [blue;red];
meanref = squeeze(nanmean(ref, 1));
stdref = squeeze(nanstd(ref, 1));
figure('Position',[0,0,800,400])
for i = 1:2
    plot(1:4, meanref(i,:),'-o', 'LineWidth', 4,'Color',col(i,:));
    hold on
end
xticks([1:4]);
xticklabels({'20%','40%','60%','80%'});
legend({'common average','bipolar'},'AutoUpdate','off')
for i = 1:2
    errorbar(1:4, meanref(i,:), stdref(i,:),'.','Color',col(i,:), 'LineWidth', 2);
end
title('Robustness of Re-reference Methods','FontSize',16,'FontWeight','bold')
ylabel('Reliability')
xlabel('% Electrode Removed')
set(gca,'FontName','Avenir','FontSize',14)
set(gca,'XLim',[0.5,4.5])
exportgraphics(gcf, strcat('figs/robustRef.png'), 'Resolution', 300);
saveas(gcf,strcat('figs/robustRef.svg'))

groupCenters = @(nGroups,nMembers,interGroupSpace) ...
    nGroups/2+.5 : nGroups+interGroupSpace : (nGroups+interGroupSpace)*nMembers-1;

figure('Position',[0,0,800,400])
boxdata{1} = squeeze(ref(:,1,:));
boxdata{2} = squeeze(ref(:,2,:));
% b = boxplotGroup(boxdata,'primaryLabels',{'',''},'SecondaryLabels',{'20%','40%','60%','80%'}, ...
%     'Colors',col,'GroupType','betweenGroups', ...
%     'PlotStyle','traditional','BoxStyle','filled', ...
%     'MedianStyle','target','Symbol','o','Widths',0.7);
hold on
b = boxplotGroup(boxdata,'primaryLabels',{'',''}, ...
    'Colors',col,'GroupType','betweenGroups', ...
    'PlotStyle','traditional','BoxStyle','outline', ...
    'Symbol','o','Widths',0.7);
title('Robustness of Re-reference Methods','FontSize',16,'FontWeight','bold')
ylabel('Reliability')
xlabel('% Electrode Removed')
ylim([0.4,1])
set(gca,'FontName','Avenir','FontSize',14)
ticks = groupCenters(numel(boxdata), size(boxdata{1},2), 1);
set(gca,'XTick',ticks,'XTickLabels',{'20%','40%','60%','80%'})
exportgraphics(gcf, strcat('figs/robustRefbox.png'), 'Resolution', 300);
saveas(gcf,strcat('figs/robustRefbox.svg'))
%% SPlit by method
for i = 1:6
    tmp = permReliability(:,sum(numFeats(1:i-1))+1:sum(numFeats(1:i)),:);
    tmp = tmp + permReliability(:,sum(numFeats(1:i-1))+1+24:sum(numFeats(1:i))+24,:);
    method(:,i,:) = nanmean(tmp,2);
    boxdata{i} = squeeze(method(:,i,:));
end
dark_cols = {'#63804f','#de7862','#51b8bd','#367068','#a5474e','#3c5488'};
dark_cols = hex2rgb(dark_cols);
dark_cols = dark_cols([4,5,6,1,2,3],:);
dark_cols = brighten(dark_cols,0.2);

figure('Position',[0,0,800,400])
% b = boxplotGroup(boxdata,'primaryLabels',{'',''},'SecondaryLabels',{'20%','40%','60%','80%'}, ...
%     'Colors',col,'GroupType','betweenGroups', ...
%     'PlotStyle','traditional','BoxStyle','filled', ...
%     'MedianStyle','target','Symbol','o','Widths',0.7);
hold on
b = boxplotGroup(boxdata,'primaryLabels',{'','','','','',''}, ...
    'Colors',dark_cols,'GroupType','betweenGroups', ...
    'PlotStyle','traditional','BoxStyle','outline', ...
    'Symbol','o','Widths',0.7);
title('Robustness of Connectivity Methods','FontSize',16,'FontWeight','bold')
ylabel('Reliability')
xlabel('% Electrode Removed')
set(gca,'FontName','Avenir','FontSize',14)
ticks = groupCenters(numel(boxdata), size(boxdata{1},2), 1);
set(gca,'XTick',ticks,'XTickLabels',{'20%','40%','60%','80%'})
exportgraphics(gcf, strcat('figs/robustMethodbox.png'), 'Resolution', 300);
saveas(gcf,strcat('figs/robustMethodbox.svg'))

%% SPlit by freq
ini = [3,10,17];
for i = 1:7
    tmp = permReliability(:,ini+i,:);
    tmp = tmp + permReliability(:,ini+i+24,:);
    freq(:,i,:) = nanmean(tmp,2);
    boxdata{i} = squeeze(freq(:,i,:));
end
dark = hex2rgb('#a5474e');
bright = hex2rgb('#ead8da');
position = [0,1];
mycolormap = interp1(position, [bright;dark], linspace(0, 1, 7), 'pchip');
figure('Position',[0,0,800,400])
% b = boxplotGroup(boxdata,'primaryLabels',{'',''},'SecondaryLabels',{'20%','40%','60%','80%'}, ...
%     'Colors',col,'GroupType','betweenGroups', ...
%     'PlotStyle','traditional','BoxStyle','filled', ...
%     'MedianStyle','target','Symbol','o','Widths',0.7);
hold on
b = boxplotGroup(boxdata,'primaryLabels',{'','','','','','',''}, ...
    'Colors',mycolormap,'GroupType','betweenGroups', ...
    'PlotStyle','traditional','BoxStyle','outline', ...
    'Symbol','o','Widths',0.7);
title('Robustness of Different Frequency Bands','FontSize',16,'FontWeight','bold')
ylabel('Reliability')
xlabel('% Electrode Removed')
set(gca,'FontName','Avenir','FontSize',14)
ticks = groupCenters(numel(boxdata), size(boxdata{1},2), 1);
set(gca,'XTick',ticks,'XTickLabels',{'20%','40%','60%','80%'})
exportgraphics(gcf, strcat('figs/robustFreqbox.png'), 'Resolution', 300);
saveas(gcf,strcat('figs/robustFreqbox.svg'))

%% T-test
results = [];
for i = 1:4
    [h,p,ci,stats] = ttest(ref(:,1,i),ref(:,2,i));
    results = [results;p];
end
