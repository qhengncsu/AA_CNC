function [Xhat, Vhat, res_norm_hist, time] = srls_GMC_matreg(y, A, d1, d2, lambda, varargin)

% [xhat, vhat, res_norm_hist] = srls_GMC_matreg(y, X, varargin)
%
% Saddle point problem:
%
% argmin_x  argmax_v { F(x,v) =
%  1/2 ||y - A x||^2 + lam ||X||_* - gamma/2 ||A(x-v)||_2^2 - lam ||V||_* }
%
% INPUT
%   y 	    response (standardized)
%   A       design matrix , each row is the vectorized matrix predictor
%   lambda  the regularization parameter
%
% OUTPUT
%   Xhat, Vhat, res_norm_hist, intercept

params = inputParser;
params.addParameter('gamma', 0.8, @(x) isnumeric(x));
params.addParameter('splitting', 'FB', @(x) ismember(x,{'DR','FB','FBF'}));
params.addParameter('max_iter', 10000, @(x) isnumeric(x));
params.addParameter('tol_stop', 1e-5, @(x) isnumeric(x));
params.addParameter('acceleration', 'aa2', @(x) ismember(x,{'original','aa2'}));
params.addParameter('early_termination', true, @(x) islogical(x));
params.addParameter('mem_size', 10, @(x) isnumeric(x));
params.addParameter('eta', 1e-2, @(x) isnumeric(x)); % before is 10^-8
params.addParameter('printevery', 100, @(x) isnumeric(x));
params.parse(varargin{:});

% soft thresholding
soft = @(x, T) max(1 - T./abs(x), 0) .* x;

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
%n = length(y);
p = d1*d2;

%
At = A';
Aty = At*y;
normA2 = norm(A)^2;
F = @(x) A*x;
FH = @(x) At*x;
if strcmp(splitting,'FB')  
    mu = 1.99*min(1,(1-gamma)/gamma)/normA2;
elseif strcmp(splitting,'FBF')
    gamma_matrix = [1-gamma,gamma;-gamma,gamma];
    mu = 0.99/(normA2*norm(gamma_matrix));
elseif strcmp(splitting,'DR')
    mu = 0.01;
    AtA = A'*A;
    G = [(1-gamma) gamma; -gamma gamma];
    M = kron(G, AtA);
    b = [Aty; zeros(p, 1)];
    left = (eye(2*p) + mu*M);
end
fprintf('step size parameter mu = %f\n', mu);

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


Xhat = reshape(xv_lambda(1:p),[d1,d2]);
Vhat = reshape(xv_lambda((p+1):(2*p)), [d1,d2]);

function zxv = forward(xv)
    x = xv(1:p,1);
    v = xv((p+1):(2*p),1);
    zx = x - mu * ( FH(F(x + gamma*(v-x))) - Aty);
    zv = v - mu * ( gamma * FH(F(v-x)) );
    zxv = [zx;zv];
end   

function zxv = forward_DR(xv)
    zxv = left\(xv+ mu*b);
end   

function xv_next = backward(zxv)
    zx = reshape(zxv(1:p),[d1,d2]);
    zv = reshape(zxv((p+1):(2*p)), [d1,d2]);
    x = svt(zx, mu * lambda);
    v = svt(zv, mu * lambda);
    xv_next = [x(:);v(:)];
end


end
