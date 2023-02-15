X = randn(1000,1000);
beta = [ones(100,1);zeros(900,1)];
y = X*beta + randn(1000,1)*std(X*beta);
groups = cell(100,1);
for i=1:100
    groups{i} = ((i-1)*10+1):(i*10);
end
%y = randn(1000,1);
y = (y-mean(y))/std(y);
X = normc(X);

tic
xhat_matrix1 = srls_GMC_path(y, X, type="single", acceleration="aa2", screen=false, gamma=0.8);
toc
tic
xhat_matrix2 = srls_GMC_path(y, X, type="single", acceleration="aa2", screen=true, gamma=0.8);
toc
norm(xhat_matrix1-xhat_matrix2)/norm(xhat_matrix2)

tic
xhat_matrix1 = srls_GMC_path(y, X, type="grouped", screen=false, lambda_min_ratio=0.01,groups=groups);
toc
tic
%acceleration can also be nesterov
xhat_matrix2 = srls_GMC_path(y, X, type="grouped",screen=true, lambda_min_ratio=0.01 ,groups=groups);
toc
norm(xhat_matrix1-xhat_matrix2)/norm(xhat_matrix2)

