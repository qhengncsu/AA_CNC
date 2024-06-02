function [Xbest, X_cell, V_cell, BICs, time] = GMC_matreg_path(y, A, d1, d2, varargin)


params = inputParser;
params.addParameter('gamma', 0.8, @(x) isnumeric(x));
params.addParameter('splitting', 'FB', @(x) ismember(x,{'DR','FB','FBF'}));
params.addParameter('max_iter', 10000, @(x) isnumeric(x));
params.addParameter('tol_stop', 1e-5, @(x) isnumeric(x));
params.addParameter('lambda_seq', double.empty(0,1), @(x) isvector(x));
%params.addParameter('nlambda', 100, @(x) isnumeric(x));
params.addParameter('acceleration', 'aa2', @(x) ismember(x,{'original','aa2'}));
params.addParameter('early_termination', true, @(x) islogical(x));
params.addParameter('mem_size', 10, @(x) isnumeric(x));
params.addParameter('eta', 1e-8, @(x) isnumeric(x));
params.addParameter('printevery', 100, @(x) isnumeric(x));
params.parse(varargin{:});



gamma = params.Results.gamma;
splitting = params.Results.splitting;
max_iter = params.Results.max_iter;
tol_stop = params.Results.tol_stop;

lambda_seq = params.Results.lambda_seq;
%nlambda = params.Results.nlambda;

acceleration = params.Results.acceleration;
early_termination = params.Results.early_termination;
mem_size = params.Results.mem_size;
eta = params.Results.eta;
printevery = params.Results.printevery;

n = length(y);

lambda0 = 0;

% Compute the LS estimate
[Xhat0, Vhat0] = srls_GMC_matreg(y, A, d1, d2, lambda0, "gamma",gamma,...
                            'splitting',splitting,'max_iter', max_iter,'tol_stop', tol_stop, ...
                             'acceleration', acceleration, 'early_termination', early_termination,...
                             'mem_size',mem_size, 'eta', eta, 'printevery', printevery);
% compute the corresponding singular values
s0 = svd(Xhat0);

% for the path
nlambda = length(lambda_seq);
BICs = NaN(nlambda, 1);
X_cell = cell(nlambda, d1, d2);
V_cell = cell(nlambda, d1, d2);

time=0;

for i=1:nlambda
    lambda = lambda_seq(i);
    % Compute the estimate at the current lambda
    [Xhat, Vhat, ~, t] = srls_GMC_matreg(y, A, d1, d2, lambda, "gamma",gamma,...
                   'splitting',splitting,'max_iter', max_iter,'tol_stop', tol_stop, ...
                   'acceleration', acceleration, 'early_termination', early_termination,...
                    'mem_size',mem_size, 'eta', eta, 'printevery', printevery);
    s = svd(Xhat);
    BICs(i) = compute_BIC(s, s0,lambda, n, d1, d2);

    X_cell{i} = Xhat;
    V_cell{i} = Vhat;
     
    time = time+t;
end


% find minimum BIC
[~, idx] = min(BICs);
Xbest = X_cell{idx};


end
