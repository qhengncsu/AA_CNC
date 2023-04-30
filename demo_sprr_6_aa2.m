load data_mice

lambda1_range = logspace(1,3,10);
lambda2_range = logspace(1,3,10);

% Initialize a matrix to store the cross-validated R2 scores
cv_scores = zeros(length(lambda1_range), length(lambda2_range));

% Define the number of folds for cross-validation
n_folds = 10;
rng(1,'philox')
n_samples = size(X, 1);
train_ratio = 2/3;
meta_train_idx = randsample(n_samples, round(train_ratio*n_samples));
meta_test_idx = setdiff(1:n_samples, meta_train_idx);
X_meta_train = X(meta_train_idx, :);
Y_meta_train = Y(meta_train_idx, :);
X_meta_test = X(meta_test_idx, :);
Y_meta_test = Y(meta_test_idx, :);

% Initialize a cvpartition object for cross-validation
cvp = cvpartition(size(X_meta_train,1), 'KFold', n_folds);
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
            X_train = X_meta_train(train_idx,:);
            Y_train = Y_meta_train(train_idx,:);
            X_val = X_meta_train(val_idx,:);
            Y_val = Y_meta_train(val_idx,:);
            
            % Fit the model on the training set with the current lambda1-lambda2 pair
            Xhat = srls_GMC_sprr(Y_train, X_train, lambda1_range(i), lambda2_range(j), 'acceleration', 'aa2', 'gamma', 0.6);
            
            % Evaluate the model on the validation set
            Y_pred_val = X_val*Xhat;
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
csvwrite('cv_scores_6_aa2.csv',cv_scores)
Xhat = srls_GMC_sprr(Y_meta_train, X_meta_train, best_lambda1, best_lambda2, 'acceleration', 'aa2', 'gamma', 0.6);
Y_pred_test = X_meta_test*Xhat;
r2_test = 1 - sum(sum((Y_meta_test - Y_pred_test).^2)) / sum(sum((Y_meta_test - mean(mean(Y_meta_test))).^2));
disp(['test R2:' num2str(r2_test)]);