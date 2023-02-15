function [xhat, vhat, res_norm_hist] = srls_GMC_acc(y, X, lambda_ratio, varargin)

% [xhat, vhat, res_norm_hist] = srls_GMC_acc(y, X, varargin)
%
% srls_GMC_path: Sparse-Regularized Least Squares with generalized MC (GMC) penalty 
%
% Saddle point problem:
%
% argmin_x  argmax_v { F(x,v) =
%  1/2 ||y - A x||^2 + lam ||x||_1 - gamma/2 ||A(x-v)||_2^2 - lam ||v||_1 }
%
% INPUT
%   y 	    response (standardized)
%   X       design matrix (columns centered and unit length)
%   lambda_ratio the ratio of lambda to lambda_max
%
% OUTPUT
%   xhat, vhat, res_norm_hist
% Algorithm: Forward-backward, Theorem 25.8 in Bauschke and Combettes(2011)
% Acceleration: Inertia, Type-II Anderson

params = inputParser;
params.addParameter('type', 'single', @(x) ischar(x)||isstring(x));
params.addParameter('groups', {}, @(x) iscell(x));
params.addParameter('gamma', 0.8, @(x) isnumeric(x));
params.addParameter('max_iter', 10000, @(x) isnumeric(x));
params.addParameter('tol_stop', 1e-5, @(x) isnumeric(x));
params.addParameter('acceleration', 'aa2', @(x) ischar(x)||isstring(x));
params.addParameter('early_termination', true, @(x) islogical(x));
params.addParameter('mem_size', 5, @(x) isnumeric(x));
params.addParameter('eta', 1e-8, @(x) isnumeric(x));
params.parse(varargin{:});

% soft thresholding
soft = @(x, T) max(1 - T./abs(x), 0) .* x;

type = params.Results.type;
groups = params.Results.groups;
ngroup = length(groups);
gamma = params.Results.gamma;
max_iter = params.Results.max_iter;
tol_stop = params.Results.tol_stop;
acceleration = params.Results.acceleration;
early_termination = params.Results.early_termination;
mem_size = params.Results.mem_size;
eta = params.Results.eta;
params_fixed = struct();
params_fixed.max_iter = max_iter;
params_fixed.tol = tol_stop;
params_fixed.early_termination = early_termination;
params_fixed.mem_size = mem_size;
params_fixed.verbose = true;
params_fixed.eta = eta;
Xt = X';
rho = norm(X)^2;
A = @(x) X*x;
AH = @(x) Xt*x;
mu = 1.99/(rho*(1-2*gamma+2*gamma^2)/(1-gamma));
n = size(X,1);
p = size(X,2);
Xty = Xt*y;
if strcmp(type,'single')
    lambda_max = max(abs(Xty));
else
    group_lens = cellfun(@(x) size(x,2), groups);
    Ks = sqrt(group_lens);
    lambda_max = max(group_norm_vec(Xty,groups)./Ks',[],'all');
end
lambda = lambda_max*lambda_ratio;
% initialization
x0 = zeros(p,1);
v0 = zeros(p,1);
xv0 = [x0;v0];
[xv_lambda, iter, res_norm_hist] = fixed_iter(xv0,@F1,params_fixed,acceleration);
fprintf('lambda = %f solved in %d iterations\n', lambda, iter);

xhat = xv_lambda(1:p);
vhat = xv_lambda((p+1):(2*p));
function xv_next = F1(xv)
    x = xv(1:p,1);
    v = xv((p+1):(2*p),1);
    zx = x - mu * ( AH(A(x + gamma*(v-x))) - Xty);
    zv = v - mu * ( gamma * AH(A(v-x)) );
    if strcmp(type,'single')
        x = soft(zx, mu * lambda);
        v = soft(zv, mu * lambda);
    else
        x = soft_group(zx, mu * lambda, groups);
        v = soft_group(zv, mu * lambda, groups);
    end
    xv_next = [x;v];
end

function norm_vec = group_norm_vec(x,groups)
    norm_vec = zeros(ngroup,1);
    for j=1:ngroup
        norm_vec(j) = norm(x(groups{j}));
    end
end

function x_prox = soft_group(x,T,groups)
    x_prox = zeros(length(x),1);
    for j=1:length(groups)
        x_prox(groups{j}) = max(1-T*sqrt(length(groups{j}))/norm(x(groups{j})),0)*x(groups{j});
    end
end
end