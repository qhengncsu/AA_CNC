X = randn(1000,5000);
beta = [ones(100,1);zeros(4900,1)];
y = X*beta + randn(1000,1)*std(X*beta);
groups = cell(200,1);
for i=1:200
    groups{i} = ((i-1)*10+1):(i*10);
end
y = randn(1000,1);
y = (y-mean(y))/std(y);
X = normc(X);

%%
tic
[xhat1, vhat1, res_norm_hist1] = srls_GMC_acc(y, X, 0.1, 'type', 'single', 'acceleration', 'original',"gamma",0.8);
toc
tic 
[xhat2, vhat2, res_norm_hist2] = srls_GMC_acc(y, X, 0.1, 'type', "single", 'acceleration', "inertia","gamma",0.8);
toc
tic
[xhat3, vhat3, res_norm_hist3] = srls_GMC_acc(y, X, 0.1, 'type', "single", 'acceleration', "aa2","gamma",0.8);
toc
figure;
plot(log(res_norm_hist1), 'k-', 'LineWidth', 1.5)
hold on; 
plot(log(res_norm_hist2), 'b-', 'LineWidth', 1.5)
hold on;
plot(log(res_norm_hist3), 'r-',  'LineWidth',1.5)
%xlabel('$iteration$','Interpreter','latex');
xlabel('iteration');
ylabel('norm of redidual');
legend('Original', 'Inertia', 'AA-II','Location', 'best','Interpreter','latex');

%%
t
[xhat1, vhat1, res_norm_hist1] = srls_GMC_acc(y, X, 0.8, 'type', "grouped", 'acceleration', "original", "groups", groups);
[xhat2, vhat2, res_norm_hist2] = srls_GMC_acc(y, X, 0.8, 'type', "grouped", 'acceleration', "inertia", "groups", groups);
[xhat3, vhat3, res_norm_hist3] = srls_GMC_acc(y, X, 0.8, 'type', "grouped", 'acceleration', "aa2", "groups", groups);
figure;
plot(log(res_norm_hist1), 'k-', 'LineWidth', 1.5)
hold on; 
plot(log(res_norm_hist2), 'b-', 'LineWidth', 1.5)
hold on;
plot(log(res_norm_hist3), 'r-',  'LineWidth',1.5)
xlabel('iteration');
ylabel('norm of redidual');
legend('Original', 'Inertia', 'AA-II','Location', 'best','Interpreter','latex');

[xhat1, vhat1, res_norm_hist1] = srls_GMC_acc(y, X, 0.05, type="single", acceleration="original", gamma=0.95);
[xhat2, vhat2, res_norm_hist2] = srls_GMC_acc(y, X, 0.05, type="single", acceleration="inertia", gamma=0.95);
[xhat3, vhat3, res_norm_hist3] = srls_GMC_acc(y, X, 0.05, type="single", acceleration="nesterov", gamma=0.95);
[xhat4, vhat4, res_norm_hist4] = srls_GMC_acc(y, X, 0.05, type="single", acceleration="aa2", gamma=0.95);
hold on 
plot(log(res_norm_hist1),'DisplayName','original')
plot(log(res_norm_hist2),'DisplayName','inertia')
plot(log(res_norm_hist3),'DisplayName','nesterov')
plot(log(res_norm_hist4),'DisplayName','aa2')
l = legend('show','Location','northwest')
hold off
