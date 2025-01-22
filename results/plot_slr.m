load cvg_single;
t=tiledlayout(2,3, 'Padding', 'none', 'TileSpacing', 'compact'); 
set(gcf,'Position',[10 10 960 640])
nexttile
plot(log(res_norm_dr), 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('DRS (%0.0f s)',t5))
hold on
plot(log(res_norm_draa), 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('AA+DRS (%0.0f s)',t6))
xlabel('iteration');
ylabel('residual norm (log scale)');
title('GMC, $\lambda=0.1\lambda_{\max}$','Interpreter','latex','FontSize',10)
l = legend('show','Location','northeast','fontsize',10)
nexttile
plot(log(res_norm_fb), 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('FBS (%0.0f s)',t1))
hold on
plot(log(res_norm_fbaa), 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('AA+FBS (%0.0f s)',t2))
xlabel('iteration');
ylabel('residual norm (log scale)');
title('GMC, $\lambda=0.1\lambda_{\max}$','Interpreter','latex','FontSize',10)
l = legend('show','Location','northeast','fontsize',10)
nexttile
plot(log(res_norm_fbf), 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('FBFS (%0.0f s)',t3))
hold on
plot(log(res_norm_fbfaa), 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('AA+FBFS (%0.0f s)',t4))
xlabel('iteration');
ylabel('residual norm (log scale)');
title('GMC, $\lambda=0.1\lambda_{\max}$','Interpreter','latex','FontSize',10)
l = legend('show','Location','northeast','fontsize',10)
load cvg_group;
nexttile
plot(log(res_norm_dr), 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('DRS (%0.0f s)',t5))
hold on
plot(log(res_norm_draa), 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('AA+DRS (%0.0f s)',t6))
xlabel('iteration');
ylabel('residual norm (log scale)');
title('Group GMC, $\lambda=0.1\lambda_{\max}$','Interpreter','latex','FontSize',10)
l = legend('show','Location','northeast','fontsize',10)
nexttile
plot(log(res_norm_fb), 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('FBS (%0.0f s)',t1))
hold on
plot(log(res_norm_fbaa), 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('AA+FBS (%0.0f s)',t2))
xlabel('iteration');
ylabel('residual norm (log scale)');
title('Group GMC, $\lambda=0.1\lambda_{\max}$','Interpreter','latex','FontSize',10)
l = legend('show','Location','northeast','fontsize',10)
nexttile
plot(log(res_norm_fbf), 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('FBFS (%0.0f s)',t3))
hold on
plot(log(res_norm_fbfaa), 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('AA+FBFS (%0.0f s)',t4))
xlabel('iteration');
ylabel('residual norm (log scale)');
title('Group GMC, $\lambda=0.1\lambda_{\max}$','Interpreter','latex','FontSize',10)
l = legend('show','Location','northeast','fontsize',10)