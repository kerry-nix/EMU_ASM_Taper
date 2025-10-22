function plot_clustergram(data, level)
global figPath connMeasures freqBands refMethods numFeats nFreq nRef
colNames = {};
for i = 1:length(connMeasures)
    if numFeats(i) == 1
        colNames = [colNames, connMeasures(i)];
    elseif numFeats(i) == nFreq
         tmpNames = strcat(connMeasures(i),{'-'},freqBands);
        colNames = [colNames, tmpNames];
    end
end
finalNames = [];
for i = 1:nRef
    finalNames = [finalNames,strcat(refMethods{i},'-',colNames)];
end
blue = hex2rgb('#3c5488');red = hex2rgb('#a5474e'); white = [1,1,1]; 
colors = [blue;white;red];
positions = [0, 0.5, 1];
mycolormap = interp1(positions, colors, linspace(0, 1, 256), 'pchip');
t.Colormap = colormap(mycolormap);
t.ColumnLabels = finalNames;
t.RowLabels = finalNames;
t.ColumnLabelsRotate = 45;
fig1 = clustergram(data);
set(fig1,t);
addTitle(fig1,'Global Connectivity Measure Grouped');
f = plot(fig1);
set(f,'FontName','Avenir');
set(gcf,'Position',[1,1,1000,1000]);
exportgraphics(gcf, strcat(figPath,'/',level,'Group.png'), 'Resolution', 300);
saveas(gcf,strcat(figPath,'/',level,'Group.svg'))
close all force