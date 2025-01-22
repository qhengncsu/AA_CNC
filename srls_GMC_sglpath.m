function [xhat_matrix, vhat_matrix, lambda_seq, res_norm_hist] = srls_GMC_sglpath(y, X, varargin)

% [xhat_matrix, vhat_matrix] = srls_GMC_path(y, X, varargin)
%
% srls_GMC_sgl: Sparse-Regularized Least Squares with generalized MC (GMC)
% penalty which computes a solution grid for sparse group lasso
%
% Saddle point problem:
%
% argmin_x  argmax_v { F(x,v) =
%  1/2 ||y - A x||^2 + lambda1 ||x||_1 + lambda2 ||x||_1,2 - 
% gamma/2 ||A(x-v)||_2^2 - lammbda1 ||v||_1 - lambda2 ||v||_1,2 }
%
% INPUT
%   y 	    response (centered)
%   X       design matrix (columns centered and standardized)
%
% OUTPUT
%   xhat_matrix, vhat_matrix
% Algorithm: Davis Yin splitting
% Acceleration: Type-II Anderson

params = inputParser;
params.addParameter('groups', {}, @(x) iscell(x));
params.addParameter('gamma', 0.8, @(x) isnumeric(x));
params.addParameter('splitting', 'DY', @(x) ismember(x,{'DY'}));
params.addParameter('max_iter', 10000, @(x) isnumeric(x));
params.addParameter('tol_stop', 1e-5, @(x) isnumeric(x));
params.addParameter('lambda_min_ratio', 0.001, @(x) isnumeric(x));
params.addParameter('nlambda', 100, @(x) isnumeric(x));
params.addParameter('lambda_seq', double.empty(0,1), @(x) isvector(x));
params.addParameter('acceleration', 'aa2', @(x) ismember(x,{'original','aa2'}));
params.addParameter('early_termination', true, @(x) islogical(x));
params.addParameter('mem_size', 10, @(x) isnumeric(x));
params.addParameter('eta', 1e-2, @(x) isnumeric(x));
params.addParameter('D', 10, @(x) isnumeric(x));
params.addParameter('printevery', 100, @(x) isnumeric(x));
params.parse(varargin{:});

% soft thresholding
soft = @(x, T) max(1 - T./abs(x), 0) .* x;

groups = params.Results.groups;
gamma = params.Results.gamma;
splitting = params.Results.splitting;
max_iter = params.Results.max_iter;
tol_stop = params.Results.tol_stop;
lambda_min_ratio = params.Results.lambda_min_ratio;
lambda_seq = params.Results.lambda_seq;
nlambda = params.Results.nlambda;
acceleration = params.Results.acceleration;
early_termination = params.Results.early_termination;
mem_size = params.Results.mem_size;
eta = params.Results.eta;
D = params.Results.D;
printevery = params.Results.printevery;
params_fixed = struct();
params_fixed.splitting = splitting;
params_fixed.max_iter = max_iter;
params_fixed.tol = tol_stop;
params_fixed.early_termination = early_termination;
params_fixed.mem_size = mem_size;
params_fixed.verbose = false;
params_fixed.eta = eta;
params_fixed.printevery = printevery;
params_fixed.D = D;
params_fixed.xi = 1e-14;

n = size(X,1);
p = size(X,2);

Xt = X';
normA2 = norm(X)^2;
A = @(x) X*x;
AH = @(x) Xt*x;
if strcmp(splitting,'DY')
    mu = 1.99*min(1,(1-gamma)/gamma)/normA2;
end
Xty = Xt*y;
if isempty(lambda_seq)
    lambda_max = max(abs(Xty));
    lambda_seq = logspace(0,log10(lambda_min_ratio),100)*lambda_max;
else
    lambda_max = max(lambda_seq);
    lambda_min = min(lambda_seq);
    nlambda = length(lambda_seq);
end

xhat_matrix = zeros(nlambda,p);
vhat_matrix = zeros(nlambda,p);

for i = 1:nlambda
    lambda1 = lambda_seq(i);
    lambda2 = lambda1/19;
    projection = @(x) soft_group(x,mu*lambda2,groups);
    params_fixed.projection = projection;
    if i>1
        xv_current = [xhat_matrix(i-1,:),vhat_matrix(i-1,:)]';
    else
        xv_current = zeros(2*p,1);
    end
    [xv_lambda, iter,res_norm_hist] = fixed_iter(xv_current,@forward,@backward,params_fixed,acceleration);
    fprintf('lambda1 = %f, lambda2 = %f solved in %d iterations\n', lambda1, lambda2, iter);
    xhat_matrix(i,:) = xv_lambda(1:p);      
    vhat_matrix(i,:) = xv_lambda((p+1):2*p); 
end

function zxv = forward(xv)
    x = xv(1:p,1);
    v = xv((p+1):(2*p),1);
    zx = x - mu * ( AH(A(x + gamma*(v-x))) - Xty);
    zv = v - mu * ( gamma * AH(A(v-x)) );
    zxv = [zx;zv];
end

function xv_next = backward(zxv)
    zx = zxv(1:p,1);
    zv = zxv((p+1):(2*p),1);
    x = soft(zx, mu * lambda1);
    v = soft(zv, mu * lambda1);
    xv_next = [x;v];
end

function x_prox = soft_group(x,T,groups)
    x_prox = zeros(length(x),1);
    for j=1:length(groups)
        x_prox(groups{j}) = max(1-T*sqrt(length(groups{j}))/norm(x(groups{j})),0)*x(groups{j});
    end
end
end
