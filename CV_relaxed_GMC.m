function CV = CV_relaxed_GMC(y, X, varargin)

% [output] = CV_relaxed_GMC(y, X, varargin)
%
% CV_relaxed_GMC: Cross-validation for relaxed GMC penalized least squares
%
% INPUT
%   y 	    response (centered and unit length)
%   X       design matrix (columns centered and unit length)
%
% OUTPUT
%  CV.xhat_matrix  the GMC estimates at a grid of lambda values using the
%                  full data set
%  CV.vhat_matrix  the GMC estimates of v at a grid of lambda values using the
%                  full data set
%  CV.intercept    estimates of intercepts at a grid of lambda values using the
%                  full data set
%  CV.lambda_seq   the grid of lambda values used for CV
%  CV.nfolds       number of folds for CV
%   CV.alpha_seq   a vector containing a grid of values for the relaxed
%                  parameter alpga
%   CV.err         an cell array of length length(alpha_seq). Each contains the cv
%                  error matrix at a given alpha
%   CV.cve         an cell array of length length(alpha_seq). Each contains the cross validation
%                  errors (a vector) at each value of lambda 
%   CV.cvse        an cell array of length length(alpha_seq). Each contains the cross validation
%                  standard errors (a vector) at each value of lambda
%   CV.min_err     an cell array of length length(alpha_seq). Each contains
%                  the min error at the given alpha.
%   CV.min_lambda  an cell array of length length(alpha_seq). Each contains
%                  the selected lambda at the given alpha.
%   CV.bestbeta    an cell array of length length(alpha_seq). Each contains
%                  the best estimate (relaxed version at the selected lambda) at the given alpha.


params = inputParser;
% parameters for CV
params.addParameter('nfolds', 5, @(x) isnumeric(x));
params.addParameter('seed', 2023, @(x) isnumeric(x));
% parameters for the algorithm
params.addParameter('type', 'single', @(x) ischar(x)||isstring(x));
params.addParameter('groups', {}, @(x) iscell(x));
params.addParameter('gamma', 0.8, @(x) isnumeric(x));
params.addParameter('splitting', 'FB', @(x) ismember(x,{'DR','FB','FBF'}));
params.addParameter('max_iter', 10000, @(x) isnumeric(x));
params.addParameter('tol_stop', 1e-5, @(x) isnumeric(x));
params.addParameter('lambda_seq', double.empty(0,1), @(x) isvector(x));
params.addParameter('alpha_seq', 0:0.5:1, @(x) isvector(x));
params.addParameter('lambda_min_ratio', 0.01, @(x) isnumeric(x));
params.addParameter('screen_off_ratio', 0.05, @(x) isnumeric(x));
params.addParameter('nlambda', 100, @(x) isnumeric(x));
params.addParameter('screen', true, @(x) islogical(x));
params.addParameter('acceleration', 'aa2', @(x) ismember(x,{'original','inertia','aa2'}));
params.addParameter('early_termination', true, @(x) islogical(x));
params.addParameter('mem_size', 5, @(x) isnumeric(x));
params.addParameter('eta', 1e-8, @(x) isnumeric(x));
params.parse(varargin{:});


nfolds = params.Results.nfolds;
seed = params.Results.seed;
type = params.Results.type;
groups = params.Results.groups;
gamma = params.Results.gamma;
splitting = params.Results.splitting;
max_iter = params.Results.max_iter;
tol_stop = params.Results.tol_stop;
lambda_min_ratio = params.Results.lambda_min_ratio;
screen_off_ratio = params.Results.screen_off_ratio;
lambda_seq = params.Results.lambda_seq;
nlambda = params.Results.nlambda;
screen = params.Results.screen;
acceleration = params.Results.acceleration;
early_termination = params.Results.early_termination;
mem_size = params.Results.mem_size;
eta = params.Results.eta;
alpha_seq = params.Results.alpha_seq;


% data standardization
[n, p] = size(X);
% center = mean(X);
% scale = sqrt(sum((X - center).^2)/n);
% X = (X - center)./scale;
% y = y -mean(y);

% Compute lambda sequence
[xhat_matrix, vhat_matrix, intercept, lambda_seq] = srls_GMC_path(y, X, 'gamma', gamma,...
              'type', type, 'groups', groups, 'lambda_seq', lambda_seq, 'splitting',splitting, ...
               'max_iter', max_iter, 'tol_stop', tol_stop,'lambda_min_ratio', lambda_min_ratio, ...
               'screen_off_ratio', screen_off_ratio, 'nlambda', nlambda,'screen', screen, ...
                'acceleration', acceleration, 'early_termination', early_termination,...
                'mem_size', mem_size, 'eta', eta);
CV.xhat_matrix = xhat_matrix;
CV.vhat_matrix = vhat_matrix;
CV.intercept = intercept;
CV.lambda_seq = lambda_seq;
CV.nfolds = nfolds;


% start CV
k  = nfolds;
CV.err = {};
cv.err = nan(nfolds, length(lambda_seq));
for i=1:k
    % get the training and test set for this fold
    testidx = i:k:n;
    trainidx = 1:n;
    trainidx(testidx)=[];
    Xtest  = X(testidx,:);
    ytest = y(testidx);
    Xtrain = X(trainidx,:);
    ytrain = y(trainidx);
    
    % fit model to training data
    [train_xmatrix, train_vmatrix, train_intercept] = srls_GMC_path(ytrain, Xtrain, 'gamma', gamma,...
        'type', type, 'groups', groups, 'lambda_seq', lambda_seq, 'splitting',splitting, ...
        'max_iter', max_iter, 'tol_stop', tol_stop,'lambda_min_ratio', lambda_min_ratio, ...
        'screen_off_ratio', screen_off_ratio, 'nlambda', nlambda,'screen', screen, ...
        'acceleration', acceleration, 'early_termination', early_termination,...
        'mem_size',mem_size, 'eta', eta);
    
    % refit
    refit_xmatrix = nan(length(lambda_seq), p);
    for j = 1:length(lambda_seq)
        % get the submatrix
        X1 = Xtrain(:, find(train_xmatrix(j, :)~=0));
        % refit
        refit = fitlm(X1,ytrain, 'Intercept', false);
        % go back to the whole vector of coefficient
        beta_hat = zeros(1, p);
        beta_hat(find(train_xmatrix(j, :)~=0))=refit.Coefficients.Estimate;
        refit_xmatrix(j, :) = beta_hat;
    end
    
    for ii = 1:length(alpha_seq)
        % fix the alpha for relaxed GMC
        alpha = alpha_seq(ii);
        % relaxed GMC with the given alpha
        relaxed_xmatrix = alpha*train_xmatrix + (1-alpha)*refit_xmatrix;
        
        % compute fit to test data
        for j=1:length(lambda_seq)
            cv.err(i,j) = norm(ytest - Xtest*relaxed_xmatrix(j, :)');
        end
        CV.err{ii, 1} = cv.err;
    end
    
end


%======== compute the best estimate using relaxed GMC and full data set
% refit; we have computed xhat_matrix using the full data set
LS_xmatrix = nan(length(lambda_seq), p);
for j = 1:length(lambda_seq)
    % get the submatrix
    X1 = X(:, find(CV.xhat_matrix(j, :)~=0));
    % refit
    refit = fitlm(X1,y, 'Intercept', false);
    % go back to the whole vector of coefficient
    beta_hat = zeros(1, p);
    beta_hat(find(CV.xhat_matrix(j, :)~=0))=refit.Coefficients.Estimate;
    LS_xmatrix(j, :) = beta_hat;
end

%Now compute other outputs
for ii = 1:length(alpha_seq)
    % fix the alpha for relaxed GMC
    alpha = alpha_seq(ii);
    cv.err = CV.err{ii,1};
    % compute cve, cvse, lambda_min etc. at the given alpha
    cv.cve  = squeeze( mean(cv.err,1) );
    cv.cvse = squeeze( std(cv.err,0,1));
    
    % find minimum error parameters
    [cv.min_err, cv.min_idx] = min(cv.cve);
    
    cv.min_lambda = lambda_seq(cv.min_idx);
    
    % the best estimate
    relaxed_xmatrix = alpha*CV.xhat_matrix + (1-alpha)*LS_xmatrix;
    cv.bestbeta = relaxed_xmatrix(cv.min_idx, :);
    
    % save all outputs at the given alpha
    CV.cve{ii, 1}=cv.cve;
    CV.cvse{ii, 1}=cv.cvse;
    CV.min_err{ii, 1}=cv.min_err;
    CV.min_idx{ii, 1}=cv.min_idx;
    CV.min_lambda{ii, 1}=cv.min_lambda;
    CV.bestbeta{ii, 1}=cv.bestbeta;
    
end


CV.alpha_seq = alpha_seq;

end
