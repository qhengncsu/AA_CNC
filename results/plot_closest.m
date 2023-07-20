load closest;
t=tiledlayout(1,2, 'Padding', 'none', 'TileSpacing', 'compact'); 
set(gcf,'Position',[10 10 800 400])
nexttile
plot(log(res_norm_hist1), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Original')
hold on
plot(log(res_norm_hist2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'AA')
xlabel('iteration');
ylabel('residual norm (log scale)');
title('p=1024')
l = legend('show','Location','northeast','fontsize',10)
nexttile
plot(times_original, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Original')
hold on
plot(times_aa2, 'r-', 'LineWidth', 1.5, 'DisplayName', 'AA')
xlabel('$$\log_{2}(p)$$','Interpreter','latex');
ylabel('time (seconds)');
title('time vs dimension')