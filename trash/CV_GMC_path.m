function cv = CV_GMC_path(y, X, varargin)

% [output] = CV_GMC(y, X, varargin)
%
% CV_GMC_path: Cross-validation for GMC penalized least squares
%
% INPUT
%   y 	    response (centered and unit length)
%   X       design matrix (columns centered and unit length)
%
% OUTPUT
%   output


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


% data standardization
[n, p] = size(X);
center = mean(X);
scale = sqrt(sum((X - center).^2)/n);
X = (X - center)./scale;
y = y -mean(y);

% Compute lambda sequence
[xhat_matrix, vhat_matrix, intercept, lambda_seq] = srls_GMC_path(y, X, 'gamma', gamma,...
              'type', type, 'groups', groups, 'lambda_seq', lambda_seq, 'splitting',splitting, ...
               'max_iter', max_iter, 'tol_stop', tol_stop,'lambda_min_ratio', lambda_min_ratio, ...
               'screen_off_ratio', screen_off_ratio, 'nlambda', nlambda,'screen', screen, ...
                'acceleration', acceleration, 'early_termination', early_termination, 'mem_size',mem_size, 'eta', eta);
cv.xhat_matrix = xhat_matrix;
cv.vhat_matrix = vhat_matrix;
cv.intercept = intercept;
cv.lambda_seq = lambda_seq;
cv.nfolds = nfolds;

% start CV
k  = nfolds;
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
    train_xmatrix = srls_GMC_path(ytrain, Xtrain, 'gamma', gamma,...
              'type', type, 'groups', groups, 'lambda_seq', lambda_seq, 'splitting',splitting, ...
               'max_iter', max_iter, 'tol_stop', tol_stop,'lambda_min_ratio', lambda_min_ratio, ...
               'screen_off_ratio', screen_off_ratio, 'nlambda', nlambda,'screen', screen, ...
                'acceleration', acceleration, 'early_termination', early_termination, 'mem_size',mem_size, 'eta', eta);
    
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
    
            
    % compute fit to test data
    for j=1:length(lambda_seq)
            cv.err(i,j) = norm(ytest - Xtest*refit_xmatrix(j, :)');
    end
end

cv.cve  = squeeze( mean(cv.err,1) ); 
cv.cvse = squeeze( std(cv.err,0,1));



% find minimum error parameters
[cv.min_err, cv.min_idx] = min(cv.cve);

cv.min_lambda = lambda_seq(cv.min_idx);
%
cv.bestbeta = cv.xhat_matrix(cv.min_idx, :);

end
