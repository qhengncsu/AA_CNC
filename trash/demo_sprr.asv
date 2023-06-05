load data_mice
tic
[Xhat,~,res_norm_hist1] = srls_GMC_sprr(Y, X, 77, 129, 'acceleration', 'original', 'gamma', 0);
time_DY = toc
tic
[Xhat_aa,~,res_norm_hist2] = srls_GMC_sprr(Y, X, 77, 129, 'acceleration', 'aa2', 'gamma', 0);
time_DY_aa2 = toc
%%
Y_pred_val = X*Xhat_aa;
r2_val = 1 - sum(sum((Y - Y_pred_val).^2)) / sum(sum((Y - mean(mean(Y))).^2));
%%

lambda1_range = logspace(1,3,10);
lambda2_range = logspace(1,3,10);
xvalues = round(lambda1_range);
yvalues = round(lambda2_range);
cv_scores = csvread("cv_scores_6_aa2.csv");
t=tiledlayout(1,2, 'Padding', 'none', 'TileSpacing', 'compact'); 
set(gcf,'Position',[10 10 1200 600])
nexttile
plot(log(res_norm_hist1), 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('original (%0.0f s)',time_DY))
hold on
plot(log(res_norm_hist2), 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('aa2 (%0.0f s)',time_DY_aa2))
xlabel('iteration');
ylabel('norm of redidual (log scale)');
title('DY, lambda1=77, lambda2=129, gamma=0.6')
l = legend('show','Location','northeast')
nexttile
heatmap(xvalues,yvalues,cv_scores)
xlabel('lambda2');
ylabel('lambda1');
title('cross-validated R2, gamma=0.6')
