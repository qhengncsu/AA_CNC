function [xhat_matrix, vhat_matrix, lambda_seq, intercept] = srls_GMC_path(y, X, varargin)

% [xhat_matrix, vhat_matrix] = srls_GMC_path(y, X, varargin)
%
% srls_GMC_path: Sparse-Regularized Least Squares with generalized MC (GMC)
% penalty which computes a solution path
%
% Saddle point problem:
%
% argmin_x  argmax_v { F(x,v) =
%  1/2 ||y - A x||^2 + lam ||x||_1 - gamma/2 ||A(x-v)||_2^2 - lam ||v||_1 }
%
% INPUT
%   y 	    response (centered and unit length)
%   X       design matrix (columns centered and unit length)
%
% OUTPUT
%   xhat_matrix, vhat_matrix
% Algorithm: Forward-backward, Theorem 25.8 in Bauschke and Combettes(2011)
% Acceleration: Inertia, Type-II Anderson

params = inputParser;
params.addParameter('type', 'single', @(x) ischar(x)||isstring(x));
params.addParameter('groups', {}, @(x) iscell(x));
params.addParameter('gamma', 0.8, @(x) isnumeric(x));
params.addParameter('splitting', 'FB', @(x) ismember(x,{'FB','FBF'}));
params.addParameter('max_iter', 10000, @(x) isnumeric(x));
params.addParameter('tol_stop', 1e-5, @(x) isnumeric(x));
params.addParameter('lambda_seq', double.empty(0,1), @(x) isvector(x));
params.addParameter('lambda_min_ratio', 1e-3, @(x) isnumeric(x));
% params.addParameter('screen_off_ratio', 0.05, @(x) isnumeric(x));
params.addParameter('nlambda', 100, @(x) isnumeric(x));
% params.addParameter('screen', true, @(x) islogical(x));
params.addParameter('acceleration', 'aa2', @(x) ismember(x,{'original','aa2'}));
params.addParameter('early_termination', true, @(x) islogical(x));
params.addParameter('mem_size', 10, @(x) isnumeric(x));
params.addParameter('eta', 1e-2, @(x) isnumeric(x));
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
lambda_min_ratio = params.Results.lambda_min_ratio;
% screen_off_ratio = params.Results.screen_off_ratio;
lambda_seq = params.Results.lambda_seq;
nlambda = params.Results.nlambda;
% screen = params.Results.screen;
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
params_fixed.verbose = false;
params_fixed.eta = eta;
params_fixed.printevery = printevery;
params_fixed.D = 10;
params_fixed.xi = 1e-14;

% data standardization
n = size(X,1);
p = size(X,2);
center = mean(X);
scale = sqrt(sum((X - center).^2)/n);
X = (X - center)./scale;
ybar = mean(y);
y = y - ybar;

Xt = X';
normA2 = norm(X)^2;
A = @(x) X*x;
AH = @(x) Xt*x;
if strcmp(splitting,'FB')
    mu = 1.99*min(1,(1-gamma)/gamma)/normA2;
elseif strcmp(splitting,'FBF')
    gamma_matrix = [1-gamma,gamma;-gamma,gamma];
    mu = 0.99/(normA2*norm(gamma_matrix));
end
Xty = Xt*y;
if isempty(lambda_seq)
    if strcmp(type,'single')
        lambda_max = max(abs(Xty));
    else
        group_lens = cellfun(@(x) size(x,2), groups);
        Ks = sqrt(group_lens);
        lambda_max = max(group_norm_vec(Xty,groups)./Ks,[],'all');
    end
    lambda_seq = logspace(0,log10(lambda_min_ratio),100)*lambda_max;
else
    group_lens = cellfun(@(x) size(x,2), groups);
    Ks = sqrt(group_lens);
    lambda_max = max(lambda_seq);
    lambda_min = min(lambda_seq);
end


% initialization
xhat_matrix = zeros(nlambda,p);
vhat_matrix = zeros(nlambda,p);
intercept = zeros(nlambda, 1);
% d_prev = zeros(p,1);
% c_prev = -Xty;

for i = 2:nlambda
    lambda = lambda_seq(i);
    xv_current = [xhat_matrix(i-1,:),vhat_matrix(i-1,:)]';
    [xv_lambda, iter] = fixed_iter(xv_current,@forward2,@backward2,params_fixed,acceleration);
    fprintf('lambda = %f solved in %d iterations\n', lambda, iter);
    % unstandardization
    bb = xv_lambda(1:p);
    xhat_matrix(i,:) = bb./scale';
    intercept(i) = ybar - center*bb;       
    vhat_matrix(i,:) = xv_lambda((p+1):2*p); 
end

function zxv = forward2(xv)
    x = xv(1:p,1);
    v = xv((p+1):(2*p),1);
    zx = x - mu * ( AH(A(x + gamma*(v-x))) - Xty);
    zv = v - mu * ( gamma * AH(A(v-x)) );
    zxv = [zx;zv];
end

function xv_next = backward2(zxv)
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
