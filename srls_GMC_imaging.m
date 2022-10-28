function [xhat, vhat, res_norm_hist] = srls_GMC_imaging(y, app, lambda, varargin)

% [xhat, vhat, res_norm_hist] = srls_GMC_imaging(y, X, varargin)
%
% srls_GMC: Sparse-Regularized Least Squares with generalized MC (GMC) penalty
%
% Saddle point problem:
%
% argmin_x  argmax_v { F(x,v) =
%  1/2 ||y - A x||^2 + lam ||x||_1 - lam/2 ||B(x-v)||_2^2 - lam ||v||_1 }
%
% INPUT
%   y 	    observed image
%   app     'deblurring' or 'inpainting'
%   lambda  penalty parameter
% OUTPUT
%   xhat, vhat, res_norm_hist
% Algorithm: Forward-backward, Theorem 25.8 in Bauschke and Combettes(2011)
% Acceleration: Nesterov with restart, Type-II Anderson
params = inputParser;
params.addParameter('H', fspecial('average',5), @(x) isnumeric(x));
params.addParameter('gamma', 0.8, @(x) isnumeric(x));
params.addParameter('max_iter', 10000, @(x) isnumeric(x));
params.addParameter('tol_stop', 1e-5, @(x) isnumeric(x));
params.addParameter('early_termination', true, @(x) islogical(x));
params.addParameter('acceleration', 'aa2', @(x) ischar(x)||isstring(x));
params.addParameter('mask', ones(size(y)), @(x) isnumeric(x));
params.parse(varargin{:});
H = params.Results.H;
gamma = params.Results.gamma;
max_iter = params.Results.max_iter;
tol_stop = params.Results.tol_stop;
early_termination = params.Results.early_termination;
acceleration = params.Results.acceleration;
mask = params.Results.mask;
params_fixed = struct();
params_fixed.max_iter = max_iter;
params_fixed.tol = tol_stop;
params_fixed.early_termination = early_termination;
params_fixed.mem_size = 5;
params_fixed.verbose = true;
n1 = size(y,1);
n2 = size(y,2);
if strcmp(app,'deblurring')
    %rho = max(abs(fft2(H,n1,n2))^2,[],'all');
    rho = 1;
    A = @(x) imfilter(x,H,'circular');
    AH = @(x) imfilter(x,H,'circular');
    AHy = AH(y);
    AHA = @(x) AH(A(x));
    B = @(x) sqrt(gamma/lambda)*imfilter(x,H,'circular');
    BH = @(x) sqrt(gamma/lambda)*imfilter(x,H,'circular');
    BHB = @(x) BH(B(x));
elseif strcmp(app,'inpainting')
    H = fspecial('average',3);
    rho = 1;
    AHA = @(x) mask.*x;
    AHy = AHA(y);
    %B = @(x) x - imfilter(x,H,'circular');
    %BHB = @(x) gamma/lambda*mask.*B(B(x));
    BHB = @(x) gamma/lambda*AHA(x);
end
xv0 = [y(:);y(:)];
if gamma>0
    mu = 1.99 / ( rho * max( 1,  gamma / (1-gamma) ) );
    [xv_lambda, iter, res_norm_hist] = fixed_iter(xv0,@F1,params_fixed,acceleration);
    xhat = reshape(xv_lambda(1:(n1*n2)),[n1,n2]);
    vhat = reshape(xv_lambda((n1*n2+1):(2*n1*n2)), [n1,n2]);
else
    mu = 1 / rho;
    [x_lambda, iter, res_norm_hist] = fixed_iter(y(:),@F2,params_fixed,acceleration);
    xhat = reshape(x_lambda,[n1,n2]);
    vhat = zeros(size(y));
end
fprintf('lambda = %f solved in %d iterations\n', lambda, iter);

function xv_next = F1(xv)
    x = reshape(xv(1:(n1*n2)),[n1,n2]);
    v = reshape(xv((n1*n2+1):(2*n1*n2)), [n1,n2]);
    zx = x - mu * (AHA(x)-AHy+lambda*BHB(v-x)) ;
    zv = v - mu * (lambda*BHB(v-x));
    x = chambolle_prox_TV_stop(zx, lambda = mu * lambda, tol=1e-3);
    v = chambolle_prox_TV_stop(zv, lambda = mu * lambda, tol=1e-3);
    xv_next = [x(:);v(:)];
end

function x_next = F2(x)
    x = reshape(x,[n1,n2]);
    zx = x - mu * (AHA(x)-AHy);
    x = chambolle_prox_TV_stop(zx, lambda = mu * lambda, tol=1e-3);
    x_next = [x(:)];
end
end