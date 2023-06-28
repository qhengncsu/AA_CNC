load biclustering;
t=tiledlayout(1,2, 'Padding', 'none', 'TileSpacing', 'compact');
set(gcf,'Position',[10 10 800 400])
nexttile
plot(log(res_norm_hist1), 'b-', 'LineWidth', 1, 'DisplayName', 'Original')
xlim([1,1700])
hold on
plot(log(res_norm_hist2), 'r-', 'LineWidth', 1, 'DisplayName', 'AA')
xlabel('iteration');
ylabel('residual norm (log scale)');
l = legend('show','Location','northeast','fontsize',10)
title('$$\lambda=10^{4.5}$$','Interpreter','latex')
nexttile
plot(4:0.5:6, times_original, 'b-', 'LineWidth', 1, 'DisplayName', 'Original')
hold on
plot(4:0.5:6, times_aa2, 'r-', 'LineWidth', 1, 'DisplayName', 'AA')
xlabel('$$\log_{10}(\lambda)$$','Interpreter','latex');
ylabel('time (seconds)');
title('time vs dimension')
