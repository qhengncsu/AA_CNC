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

t=tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'compact'); 
set(gcf,'Position',[10 10 960 960])
nexttile
plot(log(res_norm_hist_convex_original), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Original')
hold on
plot(log(res_norm_hist_convex_aa2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'AA')
xlabel('iteration');
ylabel('residual norm (log scale)');
nexttile
plot(log(res_norm_hist_cnc_original), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Original')
hold on
plot(log(res_norm_hist_cnc_aa2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'AA')
xlabel('iteration');
ylabel('residual norm (log scale)');
