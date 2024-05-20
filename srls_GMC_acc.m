function [xhat, vhat, res_norm_hist, time] = srls_GMC_acc(y, X, lambda_ratio, varargin)

% [xhat, vhat, res_norm_hist] = srls_GMC_acc(y, X, varargin)
%
% srls_GMC_path: Sparse-Regularized Least Squares with generalized MC (GMC) penalty 
%
% Saddle point problem:
%
% argmin_x  argmax_v { F(x,v) =
%  1/2 ||y - A x||^2 + lam ||x||_1 - gamma/2 ||A(x-v)||_2^2 - lam ||v||_1 }
% the ||x||_1 above may be replaced with a grouped version
%
% INPUT
%   y 	    response (standardized)
%   X       design matrix (columns centered and unit length)
%   lambda_ratio the ratio of lambda to lambda_max
%
% OUTPUT
%   xhat, vhat, res_norm_hist, intercept

params = inputParser;
params.addParameter('type', 'single', @(x) ischar(x)||isstring(x));
params.addParameter('groups', {}, @(x) iscell(x));
params.addParameter('gamma', 0.8, @(x) isnumeric(x));
params.addParameter('splitting', 'FB', @(x) ismember(x,{'DR','FB','FBF'}));
params.addParameter('max_iter', 10000, @(x) isnumeric(x));
params.addParameter('tol_stop', 1e-5, @(x) isnumeric(x));
params.addParameter('acceleration', 'aa2', @(x) ismember(x,{'original','aa2'}));
params.addParameter('early_termination', true, @(x) islogical(x));
params.addParameter('mem_size', 10, @(x) isnumeric(x));
params.addParameter('eta', 1e-8, @(x) isnumeric(x));
params.addParameter('printevery', 100, @(x) isnumeric(x));
params.parse(varargin{:});

% soft thresholding
soft = @(x, T) max(1 - T./abs(x), 0) .* x;

type = params.Results.type;
groups = params.Results.groups;
ngroup = length(groups);
gamma = params.Results.gamma;
splitting = params.Results.splitting;
max_iter = params.Results.max_iter;
tol_stop = params.Results.tol_stop;
acceleration = params.Results.acceleration;
early_termination = params.Results.early_termination;
mem_size = params.Results.mem_size;
eta = params.Results.eta;
printevery = params.Results.printevery;
params_fixed = struct();
params_fixed.splitting = splitting;
params_fixed.max_iter = max_iter;
params_fixed.tol = tol_stop;
params_fixed.early_termination = early_termination;
params_fixed.mem_size = mem_size;
params_fixed.verbose = true;
params_fixed.eta = eta;
params_fixed.printevery = printevery;
params_fixed.D = 10;
params_fixed.xi = 1e-14;

% data standardization
n = size(X,1);
p = size(X,2);
% center = mean(X);
% scale = sqrt(sum((X - center).^2)/n); 
% X = (X - center)./scale;
% y = y - mean(y);

%
Xt = X';
Xty = Xt*y;
normA2 = norm(X)^2;
A = @(x) X*x;
AH = @(x) Xt*x;
if strcmp(splitting,'FB')  
    mu = 1.99*min(1,(1-gamma)/gamma)/normA2;
elseif strcmp(splitting,'FBF')
    gamma_matrix = [1-gamma,gamma;-gamma,gamma];
    mu = 0.99/(normA2*norm(gamma_matrix));
elseif strcmp(splitting,'DR')
    mu = 0.01;
    AtA = X'*X;
    G = [(1-gamma) gamma; -gamma gamma];
    M = kron(G, AtA);
    b = [Xty; zeros(p, 1)];
    left = (eye(2*p) + mu*M);
end
fprintf('step size parameter mu = %f\n', mu);
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
tic
if  strcmp(splitting,'FB') || strcmp(splitting,'FBF') 
    [xv_lambda, iter, res_norm_hist] = fixed_iter(xv0,@forward,@backward,params_fixed,acceleration);
elseif  strcmp(splitting,'DR')
    [xv_lambda, iter, res_norm_hist] = fixed_iter(xv0,@forward_DR,@backward,params_fixed,acceleration);
end
time = toc;

fprintf('lambda = %f solved in %d iterations\n', lambda, iter);

%unstandardize the estimates 
% bb = xv_lambda(1:p);
% xhat = bb./scale';
% intercept = mean(y) - center*bb;
% vhat = xv_lambda((p+1):(2*p));

xhat = xv_lambda(1:p);
vhat = xv_lambda((p+1):(2*p));

function zxv = forward(xv)
    x = xv(1:p,1);
    v = xv((p+1):(2*p),1);
    zx = x - mu * ( AH(A(x + gamma*(v-x))) - Xty);
    zv = v - mu * ( gamma * AH(A(v-x)) );
    zxv = [zx;zv];
end   

function zxv = forward_DR(xv)
    zxv = left\(xv+ mu*b);
end   

function xv_next = backward(zxv)
    zx = zxv(1:p,1);
    zv = zxv((p+1):(2*p),1);
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
