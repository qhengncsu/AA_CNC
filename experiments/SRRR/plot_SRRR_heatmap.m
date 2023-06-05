lambda1_range = logspace(1,3,10);
lambda2_range = logspace(1,3,10);
xvalues = round(lambda1_range);
yvalues = round(lambda2_range);
cv_scores = csvread("cv_scores_DYAA.csv");
%cv_scores = csvread("cv_scores_DYAA_whole.csv");
t=tiledlayout(1,1, 'Padding', 'none', 'TileSpacing', 'compact'); 
set(gcf,'Position',[10 10 800 600])
% nexttile
figure(1);
h = heatmap(xvalues,yvalues,cv_scores);
h.XLabel = '$\lambda_2$';
h.YLabel = '$\lambda_1$';
h.title(['Cross-validated ', '$R^2$'])
h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
% xlabel('$\lambda_2$', 'Interpreter','latex')
% ylabel('lambda1');
% title('Cross-validated  R^2')