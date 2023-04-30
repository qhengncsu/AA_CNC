X = randn(1000,500);
beta = [ones(5,1);-ones(5,1);zeros(490,1)];
y = X*beta + randn(1000,1)*std(X*beta);
groups = cell(100,1);
for i=1:100
    groups{i} = ((i-1)*10+1):(i*10);
end
%y = randn(1000,1);
%y = (y-mean(y))/std(y);
%X = normc(X);

mses = zeros(50,1);
for i=1:50
    X = randn(1000,500);
    beta = [ones(5,1);-ones(5,1);zeros(490,1)];
    y = X*beta + randn(1000,1)*std(X*beta);
    [xhat_matrix1,~,intercept] = srls_GMC_path(y, X, type="single", acceleration="aa2", screen=true, gamma=0.8, splitting='FB');
    [sse, j] = min(sum((xhat_matrix1 - beta').^2,2));
    betahat = xhat_matrix1(j,:)';
    mses(i) = norm(X*betahat - X*beta)^2/1000;
end

tic
xhat_matrix2 = srls_GMC_path(y, X, type="single", acceleration="aa2", screen=false, gamma=0.9, splitting='FB');
toc
diff2 = sum((xhat_matrix2 - beta').^2,2);
norm(xhat_matrix1-xhat_matrix2)/norm(xhat_matrix2)

tic
xhat_matrix1 = srls_GMC_path(y, X, type="grouped", screen=true, lambda_min_ratio=0.01,groups=groups, splitting='FB');
toc
tic
xhat_matrix2 = srls_GMC_path(y, X, type="grouped",screen=true, lambda_min_ratio=0.01 ,groups=groups, splitting='FBF');
toc
norm(xhat_matrix1-xhat_matrix2)/norm(xhat_matrix2)



