d1 = 256;
d2 = d1;
X = zeros(d1);
row1 = 105;
row2 = 148;
X(row1:row2,40:216) = 1;
X(40:216,row1:row2) = 1;
imagesc(X);
colorbar;


X2 = checkerboard(d1/8);
imagesc(X2);
colorbar;

std_X = std(X(:));
X_noisy = X + std_X * randn(d1,d2);
mask_X = binornd(1, 0.2, size(X,1), size(X,2));
lambdas = logspace(3,-1,20);
errors_X_cnc = zeros(20,1);
tic
for i=1:20
    [xhat, vhat, res_norm_hist] = srls_GMC_mc(X_noisy.*mask_X,lambdas(i),'mask',mask_X,'gamma',0.8,'acceleration','original','splitting','FB');
    errors_X_cnc(i) = norm(xhat-X,'fro');
end
t1 = toc
tic
for i=1:20
    [xhat, vhat, res_norm_hist] = srls_GMC_mc(X_noisy.*mask_X,lambdas(i),'mask',mask_X,'gamma',0.8,'acceleration','aa2','splitting','FB');
    errors_X_cnc(i) = norm(xhat-X,'fro');
end
t2 = toc
tic
for i=1:20
    [xhat, vhat, res_norm_hist] = srls_GMC_mc(X_noisy.*mask_X,lambdas(i),'mask',mask_X,'gamma',0.8,'acceleration','original','splitting','FBF');
    errors_X_cnc(i) = norm(xhat-X,'fro');
end
t3 = toc
tic
for i=1:20
    [xhat, vhat, res_norm_hist] = srls_GMC_mc(X_noisy.*mask_X,lambdas(i),'mask',mask_X,'gamma',0.8,'acceleration','aa2','splitting','FBF');
    errors_X_cnc(i) = norm(xhat-X,'fro');
end
t4 = toc

errors_X_convex = zeros(20,1);
for i=1:20
    [xhat, vhat, res_norm_hist] = srls_GMC_mc(X_noisy.*mask_X,lambdas(i),'mask',mask_X,'gamma',0,'acceleration','aa2','splitting','FB');
    errors_X_convex(i) = norm(xhat-X,'fro');
end

std_X2 = std(X2(:));
X2_noisy = X2 + std_X2 * randn(d1,d2);
mask_X2 = binornd(1, 0.2, size(X,1), size(X,2));
lambdas = logspace(3,-1,20);
errors_X2_cnc = zeros(20,1);

tic
for i=1:20
    [xhat, vhat, res_norm_hist] = srls_GMC_mc(X2_noisy.*mask_X2,lambdas(i),'mask',mask_X2,'gamma',0.8,'acceleration','original','splitting','FB');
    errors_X2_cnc(i) = norm(xhat-X2,'fro');
end
t5 = toc

tic
for i=1:20
    [xhat, vhat, res_norm_hist] = srls_GMC_mc(X2_noisy.*mask_X2,lambdas(i),'mask',mask_X2,'gamma',0.8,'acceleration','aa2','splitting','FB');
    errors_X2_cnc(i) = norm(xhat-X2,'fro');
end
t6 = toc

tic
for i=1:20
    [xhat, vhat, res_norm_hist] = srls_GMC_mc(X2_noisy.*mask_X2,lambdas(i),'mask',mask_X2,'gamma',0.8,'acceleration','original','splitting','FBF');
    errors_X2_cnc(i) = norm(xhat-X2,'fro');
end
t7 = toc

tic
for i=1:20
    [xhat, vhat, res_norm_hist] = srls_GMC_mc(X2_noisy.*mask_X2,lambdas(i),'mask',mask_X2,'gamma',0.8,'acceleration','aa2','splitting','FBF');
    errors_X2_cnc(i) = norm(xhat-X2,'fro');
end
t8 = toc

errors_X2_convex = zeros(20,1);
for i=1:20
    [xhat, vhat, res_norm_hist] = srls_GMC_mc(X2_noisy.*mask_X2,lambdas(i),'mask',mask_X2,'gamma',0,'acceleration','aa2','splitting','FB');
    errors_X2_convex(i) = norm(xhat-X2,'fro');
end