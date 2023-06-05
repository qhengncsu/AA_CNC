load SRRR.mat
%%
t=tiledlayout(1,2, 'Padding', 'none', 'TileSpacing', 'compact'); 
set(gcf,'Position',[10 10 800 400])
nexttile
plot(log(res_norm_hist_DR), 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('Original (%0.1f s)',t_DR))
hold on
plot(log(res_norm_hist_DR_aa), 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('AA (%0.1f s)',t_DR_aa))
xlabel('iteration');
ylabel('residual norm (log scale)');
title(['DR Splitting, $\lambda_1 = $', num2str(lambda1), ', $\lambda_2 = $', num2str(lambda2)], 'Interpreter','latex')
l = legend('show','Location','northeast','fontsize',10)
nexttile
plot(log(res_norm_hist_DY), 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('Original (%0.1f s)',t_DY))
hold on
plot(log(res_norm_hist_DY_aa), 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('AA (%0.1f s)',t_DY_aa))
xlabel('iteration');
ylabel('residual norm (log scale)');
title(['DY Splitting, $\lambda_1 = $', num2str(lambda1), ', $\lambda_2 = $', num2str(lambda2)], 'Interpreter','latex')
l = legend('show','Location','northeast','fontsize',10)