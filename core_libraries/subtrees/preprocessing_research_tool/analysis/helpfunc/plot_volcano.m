function plot_volcano(pc, pval, names)
% calculate negative log10 of p-values
thres = 0.05;
logThres = -log10(thres);
logPval = -log10(pval);
sigJudge = pval<thres;
nsig = length(find(sigJudge));
blue = hex2rgb('#3c5488');red = hex2rgb('#a5474e'); 
col{1} = blue;
col{2} = red;
cols = cell2mat(col(sigJudge+1)');
% % find significantly differentially expressed genes
% sigGenes = pval < sigLevel;

% plot volcano plot
figure('Position',[0,0,700,500])
hold on
xlim([-1,1]);
ylim([0,3]);
line([-1, 1],[logThres,logThres],'Color','k','LineStyle','--','LineWidth',1)
line([0,0],[0,3],'Color','k','LineStyle','--','LineWidth',1)
scatter(pc, logPval, 400, cols, 'filled');
xlabel('Percent Change');
ylabel('-Log_{10} P');
set(gca,'FontName','Avenir','FontSize',14)
a = scatter(nan,nan,400,col{1},'filled');
b = scatter(nan,nan,400,col{2},'filled');
alpha(0.5)
legend([a,b],{'Non-Sig','Sig'},'FontSize',15);
stop = max(logPval)+0.15;
space = 0.1;
ypos = stop+space*(nsig-1):-space:stop;
ypos = ypos/3;
% label points with gene names
pcSig = pc(sigJudge)*0.5+0.5;
pvalSig = logPval(sigJudge)/3;
namesSig = names(sigJudge);
for i = 1:nsig
    annotation('textarrow',[pcSig(i)+0.03,pcSig(i)],[ypos(i),pvalSig(i)+0.03],'String',namesSig(i), ...
            'FontSize', 12, 'FontName','Avenir','HeadStyle','plain');
end
