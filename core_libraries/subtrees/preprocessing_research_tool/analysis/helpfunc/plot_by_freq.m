function [outputArg1,outputArg2] = plot_by_freq(data, level)
global freqBands figPath refNames
% color
positions = [0,1];
blue = hex2rgb('#51b8bd');red = hex2rgb('#de7862'); white = [1,1,1]; 
colors = [white;brighten(red,0.2)];
mycolormap = interp1(positions, colors, linspace(0, 1, 256), 'pchip');
h = [];
h.Colormap = colormap(mycolormap);
h.GridVisible = 'off';
h.FontSize = 12;
h.XDisplayLabels = freqBands;
h.YDisplayLabels = freqBands;
h.CellLabelColor = 'none';
close all
% plot
figure('Position',[0 0 1000 400]);
subplot(1,2,1)
fig1 = heatmap(squeeze(nanmean(data(1:size(data,1)/2,:,:),1)));
set(fig1,h);caxis([0,1])
fig1.Title = refNames{1};
set(gca,'FontName','Avenir','FontSize',14);
colorbar('off')
subplot(1,2,2)
fig1 = heatmap(squeeze(nanmean(data(1+size(data,1)/2:end,:,:),1)));
set(fig1,h);caxis([0,1])
fig1.Title = refNames{2};
set(gca,'FontName','Avenir','FontSize',14);
exportgraphics(gcf, strcat(figPath,'/freqBD',level,'.png'), 'Resolution', 300);
saveas(gcf,strcat(figPath,'/freqBD',level,'.svg'))
close all