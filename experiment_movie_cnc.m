% users = readtable('users.dat','Delimiter','::');
% movies = readtable('movies.dat','Delimiter','::');
% Y = table2array(readtable('ratings.dat','Delimiter','::'));
% n_ratings = size(Y,1);
% rng('default')
% train_idx = randsample(n_ratings, round(0.5*n_ratings));
% test_idx = setdiff(1:n_ratings, train_idx);
% Ytrain = Y(train_idx,:);
% Ytest = Y(test_idx,:);
% y = sparse(int32(Ytrain(:,1)),int32(Ytrain(:,2)),Ytrain(:,3));
% y_test = sparse(int32(Ytest(:,1)),int32(Ytest(:,2)),Ytest(:,3));
% mask = sparse(int32(Ytrain(:,1)),int32(Ytrain(:,2)),1);
% mask_test = sparse(int32(Ytest(:,1)),int32(Ytest(:,2)),1);
load movie
lambdas = 10:10:50;
mses_cnc = zeros(5,1);
for i = 1:5
    [xhat, vhat, res_norm_hist] = srls_GMC_imaging(y, 'matrix completion', lambdas(i), mask=mask, gamma=0.8, acceleration ='aa2', splitting='FBF', tol_stop=1e-4);
    mses_cnc(i) = sum(mask_test.*(xhat-y_test).^2,'all')/sum(mask_test,'all');
end
mses_cnc