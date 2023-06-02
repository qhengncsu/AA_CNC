load movie
lambdas = 40:-5:5;
mses_validation_convex = zeros(8,1);
mses_test_convex = zeros(8,1);
times_convex = zeros(8,1);
tic
[xhat, vhat, res_norm_hist] = srls_GMC_matrix(y_train, 'matrix completion', lambdas(1), 'mask', mask_train, 'gamma', 0, 'acceleration', 'original', 'splitting', 'DY', 'tol_stop',1e-4, 'lower', 1, 'upper', 5 ,'printevery', 1);
mses_validation_convex(1) = sum(mask_validation.*(xhat-y_validation).^2,'all')/sum(mask_validation,'all')
mses_test_convex(1) = sum(mask_test.*(xhat-y_test).^2,'all')/sum(mask_test,'all')
times_convex(1) = toc
for i = 2:8
    tic
    [xhat, vhat, ~] = srls_GMC_matrix(y_train, 'matrix completion', lambdas(i), 'mask', mask_train, 'gamma', 0, 'acceleration', 'original', 'splitting', 'DY', 'tol_stop',1e-4, 'lower', 1, 'upper', 5 ,'printevery', 1, 'xv0', [xhat(:);vhat(:)]);
    times_convex(i) = toc;
    mses_validation_convex(i) = sum(mask_validation.*(xhat-y_validation).^2,'all')/sum(mask_validation,'all')
    mses_test_convex(i) = sum(mask_test.*(xhat-y_test).^2,'all')/sum(mask_test,'all')
end
csvwrite('times_convex_original.csv',times_convex)
csvwrite('mses_validation_convex_original.csv',mses_validation_convex)
csvwrite('mses_test_convex_original.csv',mses_test_convex)
csvwrite('res_norm_hist_convex_original.csv',res_norm_hist)