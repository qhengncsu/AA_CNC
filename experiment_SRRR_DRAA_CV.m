%%
load("data/data_mice.mat")
%%
A = normalize(X);
Y = normalize(Y);
app = 'SRRR';
p = size(A, 2); % A is n-by-p
L = size(Y, 2); % Y is n-by-L

mu = 0; % give a random mu for SPRR, the algorithm will compute a good one
z0 = zeros(3*p*L, 1);
%%
lambda1_range = logspace(1,3,10);
lambda2_range = logspace(1,3,10);

% Initialize a matrix to store the cross-validated R2 scores
cv_scores = zeros(length(lambda1_range), length(lambda2_range));
%%
% Define the number of folds for cross-validation
n_folds = 10;
rng(1,'philox')
n_samples = size(A, 1);
train_ratio = 2/3;
meta_train_idx = randsample(n_samples, round(train_ratio*n_samples));
meta_test_idx = setdiff(1:n_samples, meta_train_idx);
A_meta_train = A(meta_train_idx, :);
Y_meta_train = Y(meta_train_idx, :);
A_meta_test = A(meta_test_idx, :);
Y_meta_test = Y(meta_test_idx, :);

%%
% Initialize a cvpartition object for cross-validation
cvp = cvpartition(size(A_meta_train,1), 'KFold', n_folds);

tic
% Perform grid search over the 2D lambda1-lambda2 grid
for i = 1:length(lambda1_range)
    for j = 1:length(lambda2_range)
        
        % Initialize a vector to store the R2 scores for the current lambda1-lambda2 pair
        r2_scores = zeros(n_folds, 1);
        
        % Iterate over the folds of cross-validation
        for k = 1:n_folds
            
            % Split the data into training and validation sets
            train_idx = cvp.training(k);
            val_idx = cvp.test(k);
            A_train = A_meta_train(train_idx,:);
            Y_train = Y_meta_train(train_idx,:);
            A_val = A_meta_train(val_idx,:);
            Y_val = Y_meta_train(val_idx,:);
            
            % Fit the model on the training set with the current lambda1-lambda2 pair
            %Xhat = srls_GMC_sprr(Y_train, X_train, lambda1_range(i), lambda2_range(j), 'acceleration', 'original', 'gamma', 0.6);
            [xhat_DR, res_norm_hist_DR] = cvx_min(A_train, Y_train(:), mu, z0, app, 'lambda1', lambda1_range(i), 'lambda2', lambda2_range(j),...
                                                 'max_iter', 1e4,  'splitting', 'DR');
            
            % Evaluate the model on the validation set
            z = ( xhat_DR(1:(p*L)) + xhat_DR((p*L+1):(2*p*L)) + xhat_DR((2*p*L+1):end))/3;
            Xhat = reshape(z, [p, L]);
            Y_pred_val = A_val*Xhat;
            r2_val = 1 - sum(sum((Y_val - Y_pred_val).^2)) / sum(sum((Y_val - mean(mean(Y_val))).^2));
            
            % Add the validation R2 score to the vector
            r2_scores(k) = r2_val;
        end
        
        % Compute the mean R2 score over the folds of cross-validation
        mean_r2 = mean(r2_scores);
        
        % Store the mean R2 score in the cv_scores matrix
        cv_scores(i,j) = mean_r2;
    end
end

% Find the lambda1 and lambda2 values with the highest cross-validated R2 score
[max_r2, max_r2_index] = max(cv_scores(:));
[max_r2_i, max_r2_j] = ind2sub(size(cv_scores), max_r2_index);
best_lambda1 = lambda1_range(max_r2_i);
best_lambda2 = lambda2_range(max_r2_j);
toc
disp(['Best lambda1: ' num2str(best_lambda1)]);
disp(['Best lambda2: ' num2str(best_lambda2)]);
disp(['Best R2 score: ' num2str(max_r2)]);
csvwrite('results/cv_scores_DRAA.csv',cv_scores)
[xhat_DR, res_norm_hist_DR] = cvx_min(A_meta_train, Y_meta_train(:), mu, z0, app, 'lambda1', best_lambda1, 'lambda2', best_lambda2,...
                                                 'max_iter', 1e4, 'splitting', 'DR');
z = ( xhat_DR(1:(p*L)) + xhat_DR((p*L+1):(2*p*L)) + xhat_DR((2*p*L+1):end))/3;
Xhat = reshape(z, [p, L]);                                       
Y_pred_test = A_meta_test*Xhat;
r2_test = 1 - sum(sum((Y_meta_test - Y_pred_test).^2)) / sum(sum((Y_meta_test - mean(mean(Y_meta_test))).^2));
disp(['test R2:' num2str(r2_test)]);

save("results/DRAA_CV.mat","best_lambda1","best_lambda2","max_r2","r2_test")