function [xhat, vhat, res_norm_hist] = srls_GMC_imaging(y, app, lambda, varargin)

% [xhat, vhat, res_norm_hist] = srls_GMC_imaging(y, app, lambda, varargin)
% INPUT
%   y 	    observed image
%   app     'deblurring' or 'inpainting' or 'matrix completion'
%   lambda  penalty parameter
% OUTPUT
%   xhat, vhat, res_norm_hist
params = inputParser;
params.addParameter('H', fspecial('average',5), @(x) isnumeric(x));
params.addParameter('gamma', 0.8, @(x) isnumeric(x));
params.addParameter('splitting', 'FBF', @(x) ismember(x,{'DR','FB','FBF'}));
params.addParameter('max_iter', 10000, @(x) isnumeric(x));
params.addParameter('tol_stop', 1e-5, @(x) isnumeric(x));
params.addParameter('early_termination', true, @(x) islogical(x));
params.addParameter('acceleration', 'aa2', @(x) ismember(x,{'original','inertia','aa2'}));
params.addParameter('mask', ones(size(y)), @(x) isnumeric(x));
params.addParameter('mem_size', 5, @(x) isnumeric(x));
params.addParameter('eta', 1e-8, @(x) isnumeric(x));
params.parse(varargin{:});
H = params.Results.H;
gamma = params.Results.gamma;
splitting = params.Results.splitting;
max_iter = params.Results.max_iter;
tol_stop = params.Results.tol_stop;
early_termination = params.Results.early_termination;
acceleration = params.Results.acceleration;
mask = params.Results.mask;
mem_size = params.Results.mem_size;
eta = params.Results.eta;
params_fixed = struct();
params_fixed.splitting = splitting;
params_fixed.max_iter = max_iter;
params_fixed.tol = tol_stop;
params_fixed.early_termination = early_termination;
params_fixed.mem_size = mem_size;
params_fixed.verbose = true;
params_fixed.eta = eta;
params_fixed.printevery = 10;
n1 = size(y,1);
n2 = size(y,2);
rho = 1;

if strcmp(app,'deblurring')
    %rho = max(abs(fft2(H,n1,n2))^2,[],'all');
    A = @(x) imfilter(x,H,'circular');
    AH = @(x) imfilter(x,H,'circular');
    AHy = AH(y);
    AHA = @(x) AH(A(x));
    %B = @(x) sqrt(gamma/lambda)*imfilter(x,H,'circular');
    %BH = @(x) sqrt(gamma/lambda)*imfilter(x,H,'circular');
    %BHB = @(x) BH(B(x));
elseif strcmp(app,'inpainting')
    AHA = @(x) mask.*x;
    AHy = AHA(y);
    %B = @(x) x - imfilter(x,H,'circular');
    %BHB = @(x) gamma/lambda*mask.*B(B(x));
    %BHB = @(x) gamma/lambda*AHA(x);
elseif strcmp(app,'matrix completion')
    AHA = @(x) mask.*x;
    AHy = AHA(y);
    %BHB = @(x) gamma/lambda*AHA(x);
end
xv0 = full([y(:);0*y(:)]);
if strcmp(splitting,'FBF') %&& strcmp(app,'matrix completion')
    gamma_matrix = [1-gamma,gamma;-gamma,gamma];
    %mu = 1.99/(rho*(1-2*gamma+2*gamma^2)/(1-gamma))
    mu = 0.99/(rho*norm(gamma_matrix))
else
    mu = 1.99/(rho*(1-2*gamma+2*gamma^2)/(1-gamma))
end
[xv_lambda, iter, res_norm_hist] = fixed_iter(xv0,@forward,@backward,params_fixed,acceleration);
xhat = reshape(xv_lambda(1:(n1*n2)),[n1,n2]);
vhat = reshape(xv_lambda((n1*n2+1):(2*n1*n2)), [n1,n2]);
fprintf('lambda = %f solved in %d iterations\n', lambda, iter);

function zxv = forward(xv)
    x = reshape(xv(1:(n1*n2)),[n1,n2]);
    v = reshape(xv((n1*n2+1):(2*n1*n2)), [n1,n2]);
    %zx = x - mu * (AHA(x)-AHy+lambda*BHB(v-x));
    zx = x - mu * (AHA(x+gamma*(v-x))-AHy);
    if gamma>0
        zv = v - mu * (gamma*AHA(v-x));
    else
        zv = v;
    end
    zxv = [zx(:);zv(:)];
end

function xv_next = backward(zxv)
    if strcmp(app,'matrix completion')
        zx = reshape(zxv(1:(n1*n2)),[n1,n2]);
        zv = reshape(zxv((n1*n2+1):(2*n1*n2)), [n1,n2]);
        x = svt(zx, mu * lambda);
        if gamma>0
        	v = svt(zv, mu * lambda);
        else
            v = zv;
        end
        xv_next = [x(:);v(:)];
    else
        zx = reshape(zxv(1:(n1*n2)),[n1,n2]);
        zv = reshape(zxv((n1*n2+1):(2*n1*n2)), [n1,n2]);
        x = chambolle_prox_TV_stop(zx, lambda = mu * lambda, maxiter = 50);
        if gamma>0
        	v = chambolle_prox_TV_stop(zv, lambda = mu * lambda, maxiter = 50);
        else
            v = zv;
        end
        xv_next = [x(:);v(:)];
    end
end
end