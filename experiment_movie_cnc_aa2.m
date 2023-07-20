load("data/movie.mat")
lambdas = 10:10:40;
mses_validation_cnc_aa2_DY = zeros(4,1);
mses_test_cnc_aa2_DY = zeros(4,1);
times_cnc_aa2_DY = zeros(4,1);
mses_validation_cnc_aa2_FB = zeros(4,1);
mses_test_cnc_aa2_FB = zeros(4,1);
times_cnc_aa2_FB = zeros(4,1);
mses_validation_cnc_aa2_FBF = zeros(4,1);
mses_test_cnc_aa2_FBF = zeros(4,1);
times_cnc_aa2_FBF = zeros(4,1);
tic
[xhat, vhat, res_norm_hist_cnc_aa2_DY] = srls_GMC_matrix(y_train, 'matrix completion', lambdas(1), 'mask', mask_train, 'gamma', 0.8, 'acceleration', 'aa2', 'splitting', 'DY', 'tol_stop',1e-4, 'lower', 1, 'upper', 5 ,'printevery', 1);
mses_validation_cnc_aa2_DY(1) = sum(mask_validation.*(xhat-y_validation).^2,'all')/sum(mask_validation,'all')
mses_test_cnc_aa2_DY(1) = sum(mask_test.*(xhat-y_test).^2,'all')/sum(mask_test,'all')
times_cnc_aa2_DY(1) = toc
for i = 2:4
    tic
    [xhat, vhat, ~] = srls_GMC_matrix(y_train, 'matrix completion', lambdas(i), 'mask', mask_train, 'gamma', 0.8, 'acceleration', 'aa2', 'splitting', 'DY', 'tol_stop',1e-4, 'lower', 1, 'upper', 5 ,'printevery', 1, 'xv0', [xhat(:);vhat(:)]);
    times_cnc_aa2_DY(i) = toc
    mses_validation_cnc_aa2_DY(i) = sum(mask_validation.*(xhat-y_validation).^2,'all')/sum(mask_validation,'all')
    mses_test_cnc_aa2_DY(i) = sum(mask_test.*(xhat-y_test).^2,'all')/sum(mask_test,'all')
end
tic
[xhat, vhat, res_norm_hist_cnc_aa2_FB] = srls_GMC_matrix(y_train, 'matrix completion', lambdas(1), 'mask', mask_train, 'gamma', 0.8, 'acceleration', 'aa2', 'splitting', 'FB', 'tol_stop',1e-4, 'lower', 1, 'upper', 5 ,'printevery', 1);
mses_validation_cnc_aa2_FB(1) = sum(mask_validation.*(xhat-y_validation).^2,'all')/sum(mask_validation,'all')
mses_test_cnc_aa2_FB(1) = sum(mask_test.*(xhat-y_test).^2,'all')/sum(mask_test,'all')
times_cnc_aa2_FB(1) = toc
for i = 2:4
    tic
    [xhat, vhat, ~] = srls_GMC_matrix(y_train, 'matrix completion', lambdas(i), 'mask', mask_train, 'gamma', 0.8, 'acceleration', 'aa2', 'splitting', 'FB', 'tol_stop',1e-4, 'lower', 1, 'upper', 5 ,'printevery', 1, 'xv0', [xhat(:);vhat(:)]);
    times_cnc_aa2_FB(i) = toc
    mses_validation_cnc_aa2_FB(i) = sum(mask_validation.*(xhat-y_validation).^2,'all')/sum(mask_validation,'all')
    mses_test_cnc_aa2_FB(i) = sum(mask_test.*(xhat-y_test).^2,'all')/sum(mask_test,'all')
end
tic
[xhat, vhat, res_norm_hist_cnc_aa2_FBF] = srls_GMC_matrix(y_train, 'matrix completion', lambdas(1), 'mask', mask_train, 'gamma', 0.8, 'acceleration', 'aa2', 'splitting', 'FBF', 'tol_stop',1e-4, 'lower', 1, 'upper', 5 ,'printevery', 1);
mses_validation_cnc_aa2_FBF(1) = sum(mask_validation.*(xhat-y_validation).^2,'all')/sum(mask_validation,'all')
mses_test_cnc_aa2_FBF(1) = sum(mask_test.*(xhat-y_test).^2,'all')/sum(mask_test,'all')
times_cnc_aa2_FBF(1) = toc
for i = 2:4
    tic
    [xhat, vhat, ~] = srls_GMC_matrix(y_train, 'matrix completion', lambdas(i), 'mask', mask_train, 'gamma', 0.8, 'acceleration', 'aa2', 'splitting', 'FBF', 'tol_stop',1e-4, 'lower', 1, 'upper', 5 ,'printevery', 1, 'xv0', [xhat(:);vhat(:)]);
    times_cnc_aa2_FBF(i) = toc
    mses_validation_cnc_aa2_FBF(i) = sum(mask_validation.*(xhat-y_validation).^2,'all')/sum(mask_validation,'all')
    mses_test_cnc_aa2_FBF(i) = sum(mask_test.*(xhat-y_test).^2,'all')/sum(mask_test,'all')
end
save("results/movie_cnc_aa2.mat","times_cnc_aa2_DY","mses_validation_cnc_aa2_DY","mses_test_cnc_aa2_DY","res_norm_hist_cnc_aa2_DY",...
    "times_cnc_aa2_FB","mses_validation_cnc_aa2_FB","mses_test_cnc_aa2_FB","res_norm_hist_cnc_aa2_FB",...
    "times_cnc_aa2_FBF","mses_validation_cnc_aa2_FBF","mses_test_cnc_aa2_FBF","res_norm_hist_cnc_aa2_FBF")