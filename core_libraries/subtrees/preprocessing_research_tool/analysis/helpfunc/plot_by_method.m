function plot_by_method(data,level)
global connMeasureNames figPath numFeats nConn refMethods
% set colors
blue = hex2rgb('#51b8bd');red = hex2rgb('#de7862'); white = [1,1,1]; 
positions = [0,0.5,1];
colors = [blue;white;red];
colors = brighten(colors,0.2);
mycolormap = interp1(positions, colors, linspace(0, 1, 256), 'pchip');
h = [];
h.Colormap = colormap(mycolormap);
h.GridVisible = 'off';
h.FontSize = 12;
h.CellLabelColor = 'none';
h.XDisplayLabels = connMeasureNames;
h.YDisplayLabels = connMeasureNames;
% data to extract
tofill = [];tmp = 0;
for i = 1:nConn
    tmp = tmp+sum(numFeats(i));
    tofill = [tofill,tmp];
end
% plot
figure('Position',[0 0 1050 400]);
subplot(1,2,1)
fig1 = heatmap(squeeze(nanmean(data(:,tofill,tofill),1)));
set(fig1,h);caxis([-1,1])
colorbar('off')
fig1.Title = refMethods{1};
set(gca,'FontName','Avenir','FontSize',14);
subplot(1,2,2)
fig1 = heatmap(squeeze(nanmean(data(:,tofill+sum(numFeats),tofill+sum(numFeats)),1)));
set(fig1,h);caxis([-1,1])
fig1.Title = refMethods{2};
set(gca,'FontName','Avenir','FontSize',14);
exportgraphics(gcf, strcat(figPath,'/methodBD',level,'.png'), 'Resolution', 300);
saveas(gcf,strcat(figPath,'/methodBD',level,'.svg'))
close all
