load("data/movie.mat")
lambdas = 5:5:15;
mses_validation_convex_original_DY = zeros(3,1);
mses_test_convex_original_DY = zeros(3,1);
times_convex_original_DY = zeros(3,1);
mses_validation_convex_original_FB = zeros(3,1);
mses_test_convex_original_FB = zeros(3,1);
times_convex_original_FB = zeros(3,1);
mses_validation_convex_original_FBF = zeros(3,1);
mses_test_convex_original_FBF = zeros(3,1);
times_convex_original_FBF = zeros(3,1);
tic
[xhat, vhat, res_norm_hist_convex_original_DY] = srls_GMC_matrix(y_train, 'matrix completion', lambdas(1), 'mask', mask_train, 'gamma', 0.8, 'acceleration', 'aa2', 'splitting', 'DY', 'tol_stop',1e-4, 'lower', 1, 'upper', 5 ,'printevery', 1);
mses_validation_convex_original_DY(1) = sum(mask_validation.*(xhat-y_validation).^2,'all')/sum(mask_validation,'all')
mses_test_convex_original_DY(1) = sum(mask_test.*(xhat-y_test).^2,'all')/sum(mask_test,'all')
times_convex_original_DY(1) = toc
for i = 2:3
    tic
    [xhat, vhat, ~] = srls_GMC_matrix(y_train, 'matrix completion', lambdas(i), 'mask', mask_train, 'gamma', 0.8, 'acceleration', 'aa2', 'splitting', 'DY', 'tol_stop',1e-4, 'lower', 1, 'upper', 5 ,'printevery', 1, 'xv0', [xhat(:);vhat(:)]);
    times_convex_original_DY(i) = toc
    mses_validation_convex_original_DY(i) = sum(mask_validation.*(xhat-y_validation).^2,'all')/sum(mask_validation,'all')
    mses_test_convex_original_DY(i) = sum(mask_test.*(xhat-y_test).^2,'all')/sum(mask_test,'all')
end
tic
[xhat, vhat, res_norm_hist_convex_original_FB] = srls_GMC_matrix(y_train, 'matrix completion', lambdas(1), 'mask', mask_train, 'gamma', 0.8, 'acceleration', 'aa2', 'splitting', 'FB', 'tol_stop',1e-4, 'lower', 1, 'upper', 5 ,'printevery', 1);
mses_validation_convex_original_FB(1) = sum(mask_validation.*(xhat-y_validation).^2,'all')/sum(mask_validation,'all')
mses_test_convex_original_FB(1) = sum(mask_test.*(xhat-y_test).^2,'all')/sum(mask_test,'all')
times_convex_original_FB(1) = toc
for i = 2:3
    tic
    [xhat, vhat, ~] = srls_GMC_matrix(y_train, 'matrix completion', lambdas(i), 'mask', mask_train, 'gamma', 0.8, 'acceleration', 'aa2', 'splitting', 'FB', 'tol_stop',1e-4, 'lower', 1, 'upper', 5 ,'printevery', 1, 'xv0', [xhat(:);vhat(:)]);
    times_convex_original_FB(i) = toc
    mses_validation_convex_original_FB(i) = sum(mask_validation.*(xhat-y_validation).^2,'all')/sum(mask_validation,'all')
    mses_test_convex_original_FB(i) = sum(mask_test.*(xhat-y_test).^2,'all')/sum(mask_test,'all')
end
tic
[xhat, vhat, res_norm_hist_convex_original_FBF] = srls_GMC_matrix(y_train, 'matrix completion', lambdas(1), 'mask', mask_train, 'gamma', 0.8, 'acceleration', 'aa2', 'splitting', 'FBF', 'tol_stop',1e-4, 'lower', 1, 'upper', 5 ,'printevery', 1);
mses_validation_convex_original_FBF(1) = sum(mask_validation.*(xhat-y_validation).^2,'all')/sum(mask_validation,'all')
mses_test_convex_original_FBF(1) = sum(mask_test.*(xhat-y_test).^2,'all')/sum(mask_test,'all')
times_convex_original_FBF(1) = toc
for i = 2:3
    tic
    [xhat, vhat, ~] = srls_GMC_matrix(y_train, 'matrix completion', lambdas(i), 'mask', mask_train, 'gamma', 0.8, 'acceleration', 'aa2', 'splitting', 'FBF', 'tol_stop',1e-4, 'lower', 1, 'upper', 5 ,'printevery', 1, 'xv0', [xhat(:);vhat(:)]);
    times_convex_original_FBF(i) = toc
    mses_validation_convex_original_FBF(i) = sum(mask_validation.*(xhat-y_validation).^2,'all')/sum(mask_validation,'all')
    mses_test_convex_original_FBF(i) = sum(mask_test.*(xhat-y_test).^2,'all')/sum(mask_test,'all')
end
save("results/movie_convex_original.mat","times_convex_original_DY","mses_validation_convex_original_DY","mses_test_convex_original_DY","res_norm_hist_convex_original_DY",...
    "times_convex_original_FB","mses_validation_convex_original_FB","mses_test_convex_original_FB","res_norm_hist_convex_original_FB",...
    "times_convex_original_FBF","mses_validation_convex_original_FBF","mses_test_convex_original_FBF","res_norm_hist_convex_original_FBF")