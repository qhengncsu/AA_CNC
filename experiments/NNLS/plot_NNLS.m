load('NNLS_results_n3e4_s5.mat')
%%
t=tiledlayout(1,2, 'Padding', 'none', 'TileSpacing', 'compact'); 
set(gcf,'Position',[10 10 800 400])
nexttile
plot(log(res_norm_hist_FB), 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('Original (%0.1f s)',t_FB))
hold on
plot(log(res_norm_hist_FB_aa), 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('AA (%0.1f s)',t_FB_aa))
xlabel('iteration');
ylabel('residual norm (log scale)');
title('FB Splitting');
l = legend('show','Location','northeast','fontsize',10)
nexttile
plot(log(res_norm_hist_DR), 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('Original (%0.1f s)',t_DR))
hold on
plot(log(res_norm_hist_DR_aa), 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('AA (%0.1f s)',t_DR_aa))
xlabel('iteration');
ylabel('residual norm (log scale)');
title('DR Splitting');
l = legend('show','Location','northeast','fontsize',10)
