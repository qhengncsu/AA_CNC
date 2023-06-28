load("data/movie.mat")
lambdas = 40:-5:5;
mses_validation_cnc_original = zeros(8,1);
mses_test_cnc_original = zeros(8,1);
times_cnc_original = zeros(8,1);
tic
[xhat, vhat, res_norm_hist_cnc_original] = srls_GMC_matrix(y_train, 'matrix completion', lambdas(1), 'mask', mask_train, 'gamma', 0.8, 'acceleration', 'original', 'splitting', 'DY', 'tol_stop',1e-4, 'lower', 1, 'upper', 5 ,'printevery', 1);
mses_validation_cnc_original(1) = sum(mask_validation.*(xhat-y_validation).^2,'all')/sum(mask_validation,'all')
mses_test_cnc_original(1) = sum(mask_test.*(xhat-y_test).^2,'all')/sum(mask_test,'all')
times_cnc_original(1) = toc
for i = 2:8
    tic
    [xhat, vhat, ~] = srls_GMC_matrix(y_train, 'matrix completion', lambdas(i), 'mask', mask_train, 'gamma', 0.8, 'acceleration', 'original', 'splitting', 'DY', 'tol_stop',1e-4, 'lower', 1, 'upper', 5 ,'printevery', 1, 'xv0', [xhat(:);vhat(:)]);
    times_cnc_original(i) = toc
    mses_validation_cnc_original(i) = sum(mask_validation.*(xhat-y_validation).^2,'all')/sum(mask_validation,'all')
    mses_test_cnc_original(i) = sum(mask_test.*(xhat-y_test).^2,'all')/sum(mask_test,'all')
end
save("results/movie_cnc_original.mat","times_cnc_original","mses_validation_cnc_original","mses_test_cnc_original","res_norm_hist_cnc_original")