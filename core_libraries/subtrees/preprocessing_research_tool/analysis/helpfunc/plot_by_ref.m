function plot_by_ref(data,error,level)
global connMeasureNames numFeats nFreq nConn figPath 
% colors
dark_cols = {'#63804f','#de7862','#51b8bd','#367068','#a5474e','#3c5488'};
bright_cols = {'#d3ddcc','#e8d0cd','#cfe2e1','#d1e0dd','#ead8da','#d0d6e2'};
dark_cols = hex2rgb(dark_cols); dark_cols = dark_cols([4,5,6,1,2,3],:);
bright_cols = hex2rgb(bright_cols); bright_cols = bright_cols([4,5,6,1,2,3],:);
bright_cols = brighten(bright_cols,0.2);
position = [0,1];
mycolormap = [];
for i = 1:length(bright_cols)
    mycolormap = [mycolormap; reshape(interp1(position, [bright_cols(i,:);dark_cols(i,:)], linspace(0, 1, 7), 'pchip'),[1,7,3])];
end
% plot
full_data = []; full_error = [];
for i = 1:nConn
    full_data = [full_data;[nan(1, nFreq-numFeats(i)) data(sum(numFeats(1:i-1))+1:sum(numFeats(1:i-1))+numFeats(i))]];
    full_error = [full_error;[nan(1, nFreq-numFeats(i)) error(sum(numFeats(1:i-1))+1:sum(numFeats(1:i-1))+numFeats(i))]];
end
bar_error_plot(full_data,full_error,mycolormap,connMeasureNames);
title(strcat(level,' Correlation between Results from Different Re-references Methods'),'Fontsize',16,'FontWeight','bold')
exportgraphics(gcf, strcat(figPath,'/refBD',level,'.png'), 'Resolution', 300);
saveas(gcf,strcat(figPath,'/refBD',level,'.svg'))
close all