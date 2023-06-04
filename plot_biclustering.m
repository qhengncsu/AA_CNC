resnorm_biclustering_aa2 = csvread('resnorm_biclustering_aa2.csv');
resnorm_biclustering_original = csvread('resnorm_biclustering_original.csv');
times_biclustering_aa2 = csvread('times_biclustering_aa2.csv');
times_biclustering_original = csvread('times_biclustering_original.csv');
t=tiledlayout(1,2, 'Padding', 'none', 'TileSpacing', 'compact');
set(gcf,'Position',[10 10 800 400])
nexttile
plot(log(resnorm_biclustering_original), 'b-', 'LineWidth', 1, 'DisplayName', 'Original')
xlim([1,1700])
hold on
plot(log(resnorm_biclustering_aa2), 'r-', 'LineWidth', 1, 'DisplayName', 'AA')
xlabel('iteration');
ylabel('residual norm (log scale)');
l = legend('show','Location','northeast','fontsize',10)
title('$$\lambda=10^{4.5}$$','Interpreter','latex')
nexttile
plot(4:0.5:6, times_biclustering_original, 'b-', 'LineWidth', 1, 'DisplayName', 'Original')
hold on
plot(4:0.5:6, times_biclustering_aa2, 'r-', 'LineWidth', 1, 'DisplayName', 'AA')
xlabel('$$\log_{10}(\lambda)$$','Interpreter','latex');
ylabel('time (seconds)');
title('time vs dimension')
