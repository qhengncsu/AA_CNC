load movie_cnc_aa2;
load movie_cnc_original;
load movie_convex_aa2;
load movie_convex_original;
lambdas = 5:5:40;
R2_val_convex_original = (1.25-mses_validation_convex_original)/1.25;
R2_val_convex_aa2 = (1.25-mses_validation_convex_aa2)/1.25;
R2_val_cnc_original = (1.25-mses_validation_cnc_original)/1.25;
R2_val_cnc_aa2 = (1.25-mses_validation_cnc_aa2)/1.25;
R2_test_convex_original = (1.25-mses_test_convex_original)/1.25;
R2_test_convex_aa2 = (1.25-mses_test_convex_aa2)/1.25;
R2_test_cnc_original = (1.25-mses_test_cnc_original)/1.25;
R2_test_cnc_aa2 = (1.25-mses_test_cnc_aa2)/1.25;
hours_cnc_original = sum(times_cnc_original)/3600;
hours_cnc_aa2 = sum(times_cnc_aa2)/3600;
hours_convex_original = sum(times_convex_original)/3600;
hours_convex_aa2 = sum(times_convex_aa2)/3600;

t=tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'compact'); 
set(gcf,'Position',[10 10 800 800])
nexttile
plot(log(res_norm_hist_convex_original), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Original')
hold on
plot(log(res_norm_hist_convex_aa2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'AA')
xlabel('iteration');
ylabel('residual norm (log scale)');
title('$$\gamma=0, \lambda=40$$','Interpreter','latex')
l = legend('show','Location','northeast','fontsize',10)
nexttile
plot(log(res_norm_hist_cnc_original), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Original')
hold on
plot(log(res_norm_hist_cnc_aa2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'AA')
xlabel('iteration');
ylabel('residual norm (log scale)');
title('$$\gamma=0.8, \lambda=40$$','Interpreter','latex')
nexttile
plot(lambdas,R2_val_convex_original(end:-1:1),'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('Convex, Original (%0.0f h)',hours_convex_original))
hold on
plot(lambdas,R2_val_convex_aa2(end:-1:1),'c-', 'LineWidth', 1.5, 'DisplayName', sprintf('Convex, AA (%0.0f h)',hours_convex_aa2))
hold on 
plot(lambdas,R2_val_cnc_original(end:-1:1),'g-', 'LineWidth', 1.5, 'DisplayName', sprintf('CNC, Original (%0.0f h)',hours_cnc_original))
hold on
plot(lambdas,R2_val_cnc_aa2(end:-1:1),'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('CNC, AA (%0.0f h)',hours_cnc_aa2))
xlabel('$$\lambda$$','Interpreter','latex');
ylabel('$$R^2$$','Interpreter','latex');
title('validation $R^2$','Interpreter','latex')
l = legend('show','Location','southwest','fontsize',10)
nexttile
plot(lambdas,R2_test_convex_original(end:-1:1),'b-', 'LineWidth', 1.5, 'DisplayName', 'Convex, Original')
hold on
plot(lambdas,R2_test_convex_aa2(end:-1:1),'c-', 'LineWidth', 1.5, 'DisplayName', 'Convex, AA')
hold on 
plot(lambdas,R2_test_cnc_original(end:-1:1),'g-', 'LineWidth', 1.5, 'DisplayName', 'CNC, Original')
hold on
plot(lambdas,R2_test_cnc_aa2(end:-1:1),'r-', 'LineWidth', 1.5, 'DisplayName', 'CNC, AA')
xlabel('$$\lambda$$','Interpreter','latex');
ylabel('$$R^2$$','Interpreter','latex');
title('test $R^2$','Interpreter','latex')
