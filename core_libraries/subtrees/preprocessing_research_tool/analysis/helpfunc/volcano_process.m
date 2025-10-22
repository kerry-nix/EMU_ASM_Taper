function volcano_process(filename)
global figPath MLPath
methods = table2cell(readtable(strcat(MLPath,'/methods.csv'),'ReadVariableNames',false));
data = table2array(readtable(strcat(MLPath,'/',filename,'.csv'),'ReadVariableNames',false));
folderName = strrep(filename,'_sozL','');
folderName = strrep(folderName,'_sozR','');
folderName = strrep(folderName,'_sozB','');
subtitle = {};
if contains(filename,'_temp')
    subtitle = [subtitle,'Temporal Lobe SOZ'];
end
if contains(filename,'_nontemp')
    subtitle = [subtitle,'Non-temporal SOZ'];
end
if contains(filename,'_good')
    subtitle = [subtitle,'Good Outcome'];
end
if contains(filename,'_sozL')
    subtitle = [subtitle,'Left SOZ'];
end
if contains(filename,'_sozR')
    subtitle = [subtitle,'Right SOZ'];
end
if contains(filename,'_sozB')
    subtitle = [subtitle,'Bilateral SOZ'];
end
tmp = strsplit(filename,'_');
if strcmp(tmp{2},'LR')
    titleText = {'Left/Right Connectivity Strength Difference',strjoin(subtitle,' - ')};
    data(:,4) = -data(:,4);
elseif strcmp(tmp{2},'SOZ')
    titleText = {'SOZ/non−SOZ Connectivity Strength Difference',strjoin(subtitle,' - ')};
end
plot_volcano(data(:,4),data(:,7),methods)
title(titleText)
exportgraphics(gcf, strcat(figPath,'/',folderName,'/',filename,'.png'), 'Resolution', 300);
saveas(gcf,strcat(figPath,'/',folderName,'/',filename,'.svg'))
close all