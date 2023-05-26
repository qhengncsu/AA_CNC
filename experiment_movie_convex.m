% users = readtable('users.dat','Delimiter','::');
% movies = readtable('movies.dat','Delimiter','::');
% Y = table2array(readtable('ratings.dat','Delimiter','::'));
% n_ratings = size(Y,1);
% rng('default')
% train_idx = randsample(n_ratings, round(0.5*n_ratings));
% test_idx = setdiff(1:n_ratings, train_idx);
% validation_idx  = randsample(test_idx,round(0.25*n_ratings));
% test_idx = setdiff(test_idx, validation_idx);
% Ytrain = Y(train_idx,:);
% Yvalidation = Y(validation_idx, :);
% Ytest = Y(test_idx,:);
% y_train = sparse(int32(Ytrain(:,1)),int32(Ytrain(:,2)),Ytrain(:,3));
% y_validation = sparse(int32(Yvalidation(:,1)),int32(Yvalidation(:,2)),Yvalidation(:,3));
% y_test = sparse(int32(Ytest(:,1)),int32(Ytest(:,2)),Ytest(:,3));
% mask_train = sparse(int32(Ytrain(:,1)),int32(Ytrain(:,2)),1);
% mask_validation = sparse(int32(Yvalidation(:,1)),int32(Yvalidation(:,2)),1);
% mask_test = sparse(int32(Ytest(:,1)),int32(Ytest(:,2)),1);
% save('movie.mat','y_train','y_validation','y_test','mask_train','mask_validation','mask_test')
load movie
sum(mask_train,'all')
lambdas = 40:-5:5;
mses_validation_convex = zeros(8,1);
mses_test_convex = zeros(8,1);
times_convex = zeros(8,1);
tic
[xhat, vhat, res_norm_hist] = srls_GMC_matrix(y_train, 'matrix completion', lambdas(1), 'mask', mask_train, 'gamma', 0, 'acceleration', 'aa2', 'splitting', 'DY', 'tol_stop',1e-4, 'lower', 1, 'upper', 5 ,'printevery', 1);
mses_validation_convex(1) = sum(mask_validation.*(xhat-y_validation).^2,'all')/sum(mask_validation,'all')
mses_test_convex(1) = sum(mask_test.*(xhat-y_test).^2,'all')/sum(mask_test,'all')
times_convex(1) = toc
for i = 2:8
    tic
    [xhat, vhat, res_norm_hist] = srls_GMC_matrix(y_train, 'matrix completion', lambdas(i), 'mask', mask_train, 'gamma', 0, 'acceleration', 'aa2', 'splitting', 'DY', 'tol_stop',1e-4, 'lower', 1, 'upper', 5 ,'printevery', 1, 'xv0', [xhat(:);vhat(:)]);
    times_convex(i) = toc
    mses_validation_convex(i) = sum(mask_validation.*(xhat-y_validation).^2,'all')/sum(mask_validation,'all')
    mses_test_convex(i) = sum(mask_test.*(xhat-y_test).^2,'all')/sum(mask_test,'all')
end
csvwrite('times_convex.csv',times_convex)
csvwrite('mses_validation_convex.csv',mses_validation_convex)
csvwrite('mses_test_convex.csv',mses_test_convex)
csvwrite('res_norm_hist_convex_aa2.csv',res_norm_hist)