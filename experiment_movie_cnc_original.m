load movie
lambdas = 40:-5:5;
mses_validation_cnc = zeros(8,1);
mses_test_cnc = zeros(8,1);
times_cnc = zeros(8,1);
tic
[xhat, vhat, res_norm_hist] = srls_GMC_matrix(y_train, 'matrix completion', lambdas(1), 'mask', mask_train, 'gamma', 0.8, 'acceleration', 'original', 'splitting', 'DY', 'tol_stop',1e-4, 'lower', 1, 'upper', 5 ,'printevery', 1);
mses_validation_cnc(1) = sum(mask_validation.*(xhat-y_validation).^2,'all')/sum(mask_validation,'all')
mses_test_cnc(1) = sum(mask_test.*(xhat-y_test).^2,'all')/sum(mask_test,'all')
times_cnc(1) = toc
for i = 2:8
    tic
    [xhat, vhat, res_norm_hist] = srls_GMC_matrix(y_train, 'matrix completion', lambdas(i), 'mask', mask_train, 'gamma', 0.8, 'acceleration', 'original', 'splitting', 'DY', 'tol_stop',1e-4, 'lower', 1, 'upper', 5 ,'printevery', 1, 'xv0', [xhat(:);vhat(:)]);
    times_cnc(i) = toc
    mses_validation_cnc(i) = sum(mask_validation.*(xhat-y_validation).^2,'all')/sum(mask_validation,'all')
    mses_test_cnc(i) = sum(mask_test.*(xhat-y_test).^2,'all')/sum(mask_test,'all')
end
csvwrite('times_cnc_original.csv',times_cnc)
csvwrite('mses_validation_cnc_original.csv',mses_validation_cnc)
csvwrite('mses_test_cnc_original.csv',mses_test_cnc)
csvwrite('res_norm_hist_cnc_original.csv',res_norm_hist)