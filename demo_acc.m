X = randn(1000,1000);
beta = [ones(100,1);zeros(900,1)];
y = X*beta + randn(1000,1)*std(X*beta);
groups = cell(100,1);
for i=1:100
    groups{i} = ((i-1)*10+1):(i*10);
end
y = randn(1000,1);
%y = (y-mean(y))/std(y);
X = normc(X);

[xhat1, vhat1, res_norm_hist1] = srls_GMC_acc(y, X, 0.8, type="single", acceleration="original");
[xhat2, vhat2, res_norm_hist2] = srls_GMC_acc(y, X, 0.8, type="single", acceleration="nesterov");
[xhat3, vhat3, res_norm_hist3] = srls_GMC_acc(y, X, 0.8, type="single", acceleration="aa2");
hold on 
plot(log(res_norm_hist1))
plot(log(res_norm_hist2))
plot(log(res_norm_hist3))
hold off

[xhat1, vhat1, res_norm_hist1] = srls_GMC_acc(y, X, 0.8, type="grouped", acceleration="original",groups=groups);
[xhat2, vhat2, res_norm_hist2] = srls_GMC_acc(y, X, 0.8, type="grouped", acceleration="nesterov",groups=groups);
[xhat3, vhat3, res_norm_hist3] = srls_GMC_acc(y, X, 0.8, type="grouped", acceleration="aa2",groups=groups);
hold on 
plot(log(res_norm_hist1))
plot(log(res_norm_hist2))
plot(log(res_norm_hist3))
hold off

