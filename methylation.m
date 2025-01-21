y_train = csvread("data/y_train.csv",1,0);
X_train = csvread("data/Z_train.csv",1,0);
n = size(X_train,1);
scale = sqrt(sum((X_train ).^2)); 
%X_train = X_train ./scale;
groups_raw = csvread("data/group.csv",1,0);
groups_unique = unique(groups_raw);
groups = cell(length(groups_unique),1);
for i=1:length(groups_unique)
    group_id = groups_unique(i);
    groups{i} = find(groups_raw==group_id);
end
lambda_seq = logspace(-0.2,-3,100)*max(abs(X_train'*y_train));

tic
xhat1_matrix = srls_GMC_sglpath(y_train, X_train, ...
    'groups',groups,'tol_stop', 1e-5,'lambda_seq',lambda_seq,'gamma',0,'acceleration','original');
t1 = toc

tic
xhat1_matrix = srls_GMC_sglpath(y_train, X_train, ...
    'groups',groups,'tol_stop', 1e-5,'lambda_seq',lambda_seq,'gamma',0,'acceleration','aa2');
t2 = toc

y_test = csvread("data/y_test.csv",1,0);
X_test = csvread("data/Z_test.csv",1,0);
beta_sd = csvread("data/Z_train_sd.csv",1,0);
R2s2 = zeros(100,1);
nonzeros1 = zeros(100,1);
for j=1:100
    beta = xhat1_matrix(j,:)'./beta_sd;
    y_pred = X_test*beta;
    R2s1(j) = corr(y_pred,y_test)^2;
    nonzeros1(j) = sum(beta~=0);
end

tic
xhat2_matrix = srls_GMC_sglpath(y_train, X_train, ...
    'groups',groups,'tol_stop', 1e-5,'lambda_seq',lambda_seq,'gamma',0.8,'acceleration','original');
t3 = toc

tic
xhat2_matrix = srls_GMC_sglpath(y_train, X_train, ...
    'groups',groups,'tol_stop', 1e-5,'lambda_seq',lambda_seq,'gamma',0.8,'acceleration','aa2');
t4 = toc

R2s2 = zeros(100,1);
nonzeros2 = zeros(100,1);
for j=1:100
    beta = xhat2_matrix(j,:)'./beta_sd;
    y_pred = X_test*beta;
    R2s2(j) = corr(y_pred,y_test)^2;
    nonzeros2(j) = sum(beta~=0);
end

R2s1(isnan(R2s1)) = 0;
R2s2(isnan(R2s2)) = 0;
ratio = logspace(-0.2,-3,100);
plot(ratio,R2s1, 'b-', 'LineWidth', 1.5, 'DisplayName', '$$\gamma = 0$$')
hold on
plot(ratio,R2s2, 'r-', 'LineWidth', 1.5, 'DisplayName', '$$\gamma = 0.8$$')
xlim([min(ratio)-0.01,max(ratio)+0.01])
set (gca, 'xdir', 'reverse')
l = legend('show','Location','east','Interpreter','latex')
xlabel('$$\lambda/\lambda_{\max}$$','Interpreter','latex');
ylabel('$$R^2$$','Interpreter','latex');

title("$$R^2$$ - Sparse Group Lasso",'Interpreter','latex')

[~,~,~,res_norm_hist1] = srls_GMC_sglpath(y_train, X_train, ...
    'groups',groups,'tol_stop', 1e-5,'lambda_seq',lambda_seq(1),'gamma',0.8,'acceleration','original');
[~,~,~,res_norm_hist2] = srls_GMC_sglpath(y_train, X_train, ...
    'groups',groups,'tol_stop', 1e-5,'lambda_seq',lambda_seq(1),'gamma',0.8,'acceleration','aa2','eta',0,'D',1e7);
[~,~,~,res_norm_hist3] = srls_GMC_sglpath(y_train, X_train, ...
    'groups',groups,'tol_stop', 1e-5,'lambda_seq',lambda_seq(1),'gamma',0.8,'acceleration','aa2');

plot(log(res_norm_hist1), 'b-', 'LineWidth', 1, 'DisplayName', 'Original')
hold on
plot(log(res_norm_hist2), 'g-', 'LineWidth', 1, 'DisplayName', 'Naive AA')
hold on
plot(log(res_norm_hist3), 'r-', 'LineWidth', 1, 'DisplayName', 'Algorithm 1')
xlabel('iteration');
ylabel('residual norm (log scale)');
title('DYS, $\lambda=10^{-0.2}\lambda_{\max}$','Interpreter','latex','FontSize',10)
l = legend('show','Location','northeast','fontsize',10)