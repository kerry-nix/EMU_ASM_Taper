function plot_corrs(data,titles,level)
global nConn connMeasureNames freqBands numFeats figPath
if exist(strcat(figPath,'/corr'),'dir') == 0
    mkdir(strcat(figPath,'/corr'))
end
% heatmap settings
colNames = {};
for i = 1:nConn
    if numFeats(i) == 1
        colNames = [colNames, connMeasureNames(i)];
    elseif numFeats(i) == 7
        colNames = [colNames, freqBands];
    end
end
blue = hex2rgb('#3c5488');red = hex2rgb('#a5474e'); white = [1,1,1]; 
colors = [blue;white;red];
positions = [0, 0.5, 1];
mycolormap = interp1(positions, colors, linspace(0, 1, 256), 'pchip');
% Customize the appearance of the heatmap
h.Colormap = colormap(mycolormap);%brighten(redbluecmap, 0.6);
h.ColorLimits = [-1,1];
h.GridVisible = 'off';
h.FontSize = 12;
h.XDisplayLabels = [colNames,colNames];
h.YDisplayLabels = [colNames,colNames];
h.Position = [0.125 0.15 0.5 0.825];
close all
figure('Position',[1,1,1440,796])
% annotation
annotation('textbox',[0.175, 0.03,0.1,0.075], 'String','Coherence', 'FontSize', 12,'LineStyle','none');
annotation('line',[0.165, 0.23],[0.10,0.10]);
annotation('textbox',[0.255, 0.03,0.1,0.075], 'String','PLV', 'FontSize', 12,'LineStyle','none');
annotation('line',[0.235, 0.3],[0.10,0.10]);
annotation('textbox',[0.32, 0.03,0.1,0.075], 'String',{'Relative','Entropy'}, 'FontSize', 12,'LineStyle','none');
annotation('line',[0.305, 0.37],[0.10,0.10]);
plus = 0.25;
annotation('textbox',[0.175+plus, 0.03,0.1,0.075], 'String','Coherence', 'FontSize', 12,'LineStyle','none');
annotation('line',[0.165+plus, 0.23+plus],[0.10,0.10]);
annotation('textbox',[0.255+plus, 0.03,0.1,0.075], 'String','PLV', 'FontSize', 12,'LineStyle','none');
annotation('line',[0.235+plus, 0.3+plus],[0.10,0.10]);
annotation('textbox',[0.32+plus, 0.03,0.1,0.075], 'String',{'Relative','Entropy'}, 'FontSize', 12,'LineStyle','none');
annotation('line',[0.305+plus, 0.37+plus],[0.10,0.10]);
annotation('textbox',[0.17, 0.03,0.9,0.02], 'String','Common Average Reference', 'FontSize', 12,'LineStyle','none');
annotation('line',[0.125 0.37],[0.05,0.05]);
annotation('textbox',[0.2+plus, 0.03,0.1,0.02], 'String','Bipolar Reference', 'FontSize', 12,'LineStyle','none');
annotation('line',[0.125+plus 0.37+plus],[0.05,0.05]);
% vertical
annotation('textbox',[0.05, 0.03,0.75,0.84], 'String','Coherence', 'FontSize', 12,'LineStyle','none');
annotation('line',[0.095, 0.095],[0.8,0.91]);
annotation('textbox',[0.075, 0.03,0.635,0.725], 'String','PLV', 'FontSize', 12,'LineStyle','none');
annotation('line',[0.095, 0.095],[0.685,0.795]);
annotation('textbox',[0.06, 0.03,0.55,0.61], 'String',{'Relative','Entropy'}, 'FontSize', 12,'LineStyle','none');
annotation('line',[0.095, 0.095],[0.565,0.68]);
minus = 0.415;
annotation('textbox',[0.05, 0.03,0.75-minus,0.84-minus], 'String','Coherence', 'FontSize', 12,'LineStyle','none');
annotation('line',[0.095, 0.095],[0.8-minus,0.91-minus]);
annotation('textbox',[0.075, 0.03,0.635-minus,0.725-minus], 'String','PLV', 'FontSize', 12,'LineStyle','none');
annotation('line',[0.095, 0.095],[0.685-minus,0.795-minus]);
annotation('textbox',[0.06, 0.03,0.55-minus,0.61-minus], 'String',{'Relative','Entropy'}, 'FontSize', 12,'LineStyle','none');
annotation('line',[0.095, 0.095],[0.565-minus,0.68-minus]);
annotation('textbox',[0.005, 0.03,0.67,0.76], 'String',{'Common','Average','Reference'}, 'FontSize', 12,'LineStyle','none');
annotation('line',[0.05, 0.05],[0.565,0.97]);
annotation('textbox',[0.005, 0.03,0.67-minus,0.76-minus], 'String',{'Bipolar','Reference'}, 'FontSize', 12,'LineStyle','none');
annotation('line',[0.05, 0.05],[0.565-minus,0.97-minus]);

% plotting
if strcmp('Global',level)
    fig1 = heatmap(data);
    set(fig1,h);
    fig1.Title = strrep(strcat('Global Connectivity Measure Correlation'),'_','-');
    exportgraphics(gcf, strcat(figPath,'/corr/Global.png'), 'Resolution', 300);
    saveas(gcf,strcat(figPath,'/corr/Global.svg'))
    close all
else
    for i = 1:length(titles)
        % Export the heatmap to a file
        % network level
        fig1 = heatmap(squeeze(data(i,:,:)));
        set(fig1,h);
        fig1.Title = strrep(strcat(titles{i}," ",level,' Connectivity Measure Correlation'),'_','-');
        exportgraphics(gcf, strcat(figPath,'/corr/',titles{i},'_',level,'.png'), 'Resolution', 300);
        saveas(gcf,strcat(figPath,'/corr/',titles{i},'_',level,'.svg'))
    end
    fig1 = heatmap(squeeze(mean(data,1,'omitnan')));
    set(fig1,h);
    fig1.Title = strrep(strcat("Averaged ",level,' Connectivity Measure Correlation'),'_','-');
    exportgraphics(gcf, strcat(figPath,'/corr/',level,'.png'), 'Resolution', 300);
    saveas(gcf,strcat(figPath,'/corr/',level,'.svg'))
    close all
end