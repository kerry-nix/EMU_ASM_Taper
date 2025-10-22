function plot_soz_scatter(nodeStrAll,numlabel,folderName)
global figPath finalNames
blue = hex2rgb('#51b8bd');red = hex2rgb('#de7862');green = hex2rgb('#63804f');
col{1} = blue';col{2} = red';col{3} = green';
cols = [blue;red;green];
cols = brighten(cols,0.2);
labels = unique(numlabel);
leg = {'Left SOZ','Right SOZ','Bilateral SOZ'};
leg = leg(labels);
tmp = strsplit(folderName,'_');
if strcmp(tmp{2},'LR')
    axisLabel = {'Left','Right'};
elseif strcmp(tmp{2},'SOZ')
    axisLabel = {'SOZ','Non-SOZ'};
end
if exist(strcat(figPath,'/',folderName),'dir') == 0
    mkdir(strcat(figPath,'/',folderName))
end
f = figure('Position',[0 0 700 500]);
for i = 1:size(nodeStrAll,1)
    for j = 1:length(labels)
        scatter(nodeStrAll(i,numlabel==j,1),nodeStrAll(i,numlabel==j,2),120,cols(j,:),'filled','MarkerFaceAlpha',.7);
        hold on
    end
    hold off
    title(finalNames{i})
    legend(leg,'Location','SouthEast','AutoUpdate','off')
    set(gca,'FontName','Avenir','FontSize',14)
    xlabel(axisLabel{1})
    ylabel(axisLabel{2})
    f.CurrentAxes.XLim(2) = max(f.CurrentAxes.XLim(2),f.CurrentAxes.YLim(2));
    f.CurrentAxes.YLim(2) = max(f.CurrentAxes.XLim(2),f.CurrentAxes.YLim(2));
    f.CurrentAxes.XLim(1) = min(f.CurrentAxes.XLim(1),f.CurrentAxes.YLim(1));
    f.CurrentAxes.YLim(1) = min(f.CurrentAxes.XLim(1),f.CurrentAxes.YLim(1));
    line(f.CurrentAxes.XLim,f.CurrentAxes.YLim,'Color','black','LineStyle','--')
    exportgraphics(gcf, strcat(figPath,'/',folderName,'/',num2str(i),'.png'), 'Resolution', 300);
    saveas(gcf,strcat(figPath,'/',folderName,'/',num2str(i),'.svg'))
end
close all