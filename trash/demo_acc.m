X = randn(1000,1000);
beta = [ones(100,1);zeros(900,1)];
y = X*beta + 0.5*randn(1000,1)*std(X*beta);
groups = cell(10,1);
for i=1:10
    groups{i} = ((i-1)*100+1):(i*100);
end
%y = randn(1000,1);
[n, p] = size(X);
center = mean(X);
scale = sqrt(sum((X - center).^2)/n);
X = (X - center)./scale;
y = y -mean(y);

%%
tic
[xhat1, vhat1, res_norm_hist1] = srls_GMC_acc(y, X, 0.1, 'type', 'single', 'acceleration', 'original',"splitting","FB",'gamma',0.8);
toc
tic 
[xhat2, vhat2, res_norm_hist2] = srls_GMC_acc(y, X, 0.1, 'type', 'single', 'acceleration', "original","splitting","FBF",'gamma',0.8);
toc
tic 
[xhat3, vhat3, res_norm_hist3] = srls_GMC_acc(y, X, 0.1, 'type', 'single', 'acceleration', "aa2","splitting","FB",'gamma',0.8);
toc
tic 
[xhat4, vhat4, res_norm_hist4] = srls_GMC_acc(y, X, 0.1, 'type', 'single', 'acceleration', "aa2","splitting","FBF",'gamma',0.8);
toc
figure;
plot(log(res_norm_hist1), 'c-', 'LineWidth', 1.5)
hold on; 
plot(log(res_norm_hist2), 'b-', 'LineWidth', 1.5)
hold on;
plot(log(res_norm_hist3), 'r-', 'LineWidth', 1.5)
hold on;
plot(log(res_norm_hist4), 'g-', 'LineWidth', 1.5)
hold off;
%xlabel('$iteration$','Interpreter','latex');
xlabel('iteration');
ylabel('norm of redidual');

%%
tic
[xhat1, vhat1, res_norm_hist1] = srls_GMC_acc(y, X, 0.1, 'type', 'grouped', 'acceleration', 'original',"splitting","FB",'groups',groups,'gamma',0.8);
toc
tic 
[xhat2, vhat2, res_norm_hist2] = srls_GMC_acc(y, X, 0.1, 'type', "grouped", 'acceleration', "original","splitting","FBF",'groups',groups,'gamma',0.8);
toc
tic 
[xhat3, vhat3, res_norm_hist3] = srls_GMC_acc(y, X, 0.1, 'type', "grouped", 'acceleration', "aa2","splitting","FB",'groups',groups,'gamma',0.8);
toc
tic 
[xhat4, vhat4, res_norm_hist4] = srls_GMC_acc(y, X, 0.1, 'type', "grouped", 'acceleration', "aa2","splitting","FBF",'groups',groups,'gamma',0.8);
toc
figure;
plot(log(res_norm_hist1), 'c-', 'LineWidth', 1.5)
hold on; 
plot(log(res_norm_hist2), 'b-', 'LineWidth', 1.5)
hold on;
plot(log(res_norm_hist3), 'r-', 'LineWidth', 1.5)
hold on;
plot(log(res_norm_hist4), 'g-', 'LineWidth', 1.5)
hold off;
%xlabel('$iteration$','Interpreter','latex');
xlabel('iteration');
ylabel('norm of redidual');

