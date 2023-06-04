res_norm_hist_cnc_original = csvread('res_norm_hist_cnc_original.csv');
res_norm_hist_cnc_aa2 = csvread('res_norm_hist_cnc_aa2.csv');
res_norm_hist_convex_original = csvread('res_norm_hist_convex_original.csv');
res_norm_hist_convex_aa2 = csvread('res_norm_hist_convex_aa2.csv');
mses_validation_cnc = csvread('mses_validation_cnc.csv');
mses_validation_cnc_original = csvread('mses_validation_cnc_original.csv');
mses_validation_convex = csvread('mses_validation_convex.csv');
mses_validation_convex_original = csvread('mses_validation_convex_original.csv');
mses_test_cnc = csvread('mses_test_cnc.csv');
mses_test_cnc_original = csvread('mses_test_cnc_original.csv');
mses_test_convex = csvread('mses_test_convex.csv');
mses_test_convex_original = csvread('mses_test_convex_original.csv');
lambdas = 5:5:40;
R2_val_convex_original = (1.25-mses_validation_convex_original)/1.25;
R2_val_convex = (1.25-mses_validation_convex)/1.25;
R2_val_cnc_original = (1.25-mses_validation_cnc_original)/1.25;
R2_val_cnc = (1.25-mses_validation_cnc)/1.25;
R2_test_convex_original = (1.25-mses_test_convex_original)/1.25;
R2_test_convex = (1.25-mses_test_convex)/1.25;
R2_test_cnc_original = (1.25-mses_test_cnc_original)/1.25;
R2_test_cnc = (1.25-mses_test_cnc)/1.25;

times_convex = csvread('times_convex.csv');
times_convex_original = csvread('times_convex_original.csv');
times_cnc = csvread('times_cnc.csv');
times_cnc_original = csvread('times_cnc_original.csv');

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
plot(lambdas,R2_val_convex_original(end:-1:1),'b-', 'LineWidth', 1.5, 'DisplayName', 'Convex, Original')
hold on
plot(lambdas,R2_val_convex(end:-1:1),'c-', 'LineWidth', 1.5, 'DisplayName', 'Convex, AA')
hold on 
plot(lambdas,R2_val_cnc_original(end:-1:1),'g-', 'LineWidth', 1.5, 'DisplayName', 'CNC, Original')
hold on
plot(lambdas,R2_val_cnc(end:-1:1),'r-', 'LineWidth', 1.5, 'DisplayName', 'CNC, AA')
xlabel('$$\lambda$$','Interpreter','latex');
ylabel('$$R^2$$','Interpreter','latex');
title('validation $R^2$','Interpreter','latex')
l = legend('show','Location','southwest','fontsize',10)
nexttile
plot(lambdas,R2_test_convex_original(end:-1:1),'b-', 'LineWidth', 1.5, 'DisplayName', 'Convex, Original')
hold on
plot(lambdas,R2_test_convex(end:-1:1),'c-', 'LineWidth', 1.5, 'DisplayName', 'Convex, AA')
hold on 
plot(lambdas,R2_test_cnc_original(end:-1:1),'g-', 'LineWidth', 1.5, 'DisplayName', 'CNC, Original')
hold on
plot(lambdas,R2_test_cnc(end:-1:1),'r-', 'LineWidth', 1.5, 'DisplayName', 'CNC, AA')
xlabel('$$\lambda$$','Interpreter','latex');
ylabel('$$R^2$$','Interpreter','latex');
title('test $R^2$','Interpreter','latex')
