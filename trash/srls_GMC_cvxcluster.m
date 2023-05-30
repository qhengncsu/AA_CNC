function [Xhat, Vhat, res_norm_hist] = srls_GMC_cvxcluster(U, lambda, varargin)
% [xhat, vhat, res_norm_hist] = srls_GMC_imaging(U, lambda, varargin)
% INPUT
%   y 	    observed image
%   app     'deblurring' or 'inpainting' or 'matrix completion'
%   lambda  penalty parameter
% OUTPUT
%   xhat, vhat, res_norm_hist
params = inputParser;
params.addParameter('gamma', 0.8, @(x) isnumeric(x));
params.addParameter('splitting', 'DY', @(x) ismember(x,{'DR','DY'}));
params.addParameter('max_iter', 10000, @(x) isnumeric(x));
params.addParameter('tol_stop', 1e-5, @(x) isnumeric(x));
params.addParameter('early_termination', true, @(x) islogical(x));
params.addParameter('acceleration', 'aa2', @(x) ismember(x,{'original','aa2'}));
params.addParameter('mem_size', 5, @(x) isnumeric(x));
params.addParameter('eta', 1e-8, @(x) isnumeric(x));
params.addParameter('K', 5, @(x) isnumeric(x));
params.addParameter('phi', 0, @(x) isnumeric(x));
params.parse(varargin{:});
gamma = params.Results.gamma;
splitting = params.Results.splitting;
max_iter = params.Results.max_iter;
tol_stop = params.Results.tol_stop;
early_termination = params.Results.early_termination;
acceleration = params.Results.acceleration;
mem_size = params.Results.mem_size;
eta = params.Results.eta;
K = params.Results.K;
phi = params.Results.phi;
params_fixed = struct();
params_fixed.splitting = splitting;
params_fixed.max_iter = max_iter;
params_fixed.tol = tol_stop;
params_fixed.early_termination = early_termination;
params_fixed.mem_size = mem_size;
params_fixed.verbose = true;
params_fixed.eta = eta;
mu = 1.99/((1-2*gamma+2*gamma^2)/(1-gamma));
[U, D, w] = preprocessU(U,K,phi);
AHy = U(:);
X0 = U;
V0 = 0.*U;
XD0 = D*X0;
VD0 = D*V0;
p = size(X0,2);
pX = size(X0,1)*size(X0,2);
pXD = size(XD0,1)*size(XD0,2);
z0 = [X0(:);V0(:);XD0(:);VD0(:)];
E = [kron(speye(p),D),sparse(pXD,pX),-speye(pXD),sparse(pXD,pXD);
     sparse(pXD,pX),kron(speye(p),D),sparse(pXD,pXD),-speye(pXD)];
EET = E*E';
function z_proj = projection(z)
    z_proj = z - E'* (EET\(E*z));
end
params_fixed.projection = @projection;
[z_lambda, iter, res_norm_hist] = fixed_iter(z0,@forward,@backward,params_fixed,acceleration);
Xhat = reshape(z_lambda(1:pX),size(X0));
Vhat = reshape(z_lambda((pX+1):2*pX),size(X0));
fprintf('lambda = %f solved in %d iterations\n', lambda, iter);

function z_forward = forward(z)
    X = z(1:pX);
    V = z((pX+1):2*pX);
    zX = mu * ((X+gamma*(V-X))-AHy);
    zV = mu * (gamma*(V-X));
    z_forward = [zX;zV;z((2*pX+1):end)];
end

function z_backward = backward(z)
    XD = reshape(z((2*pX+1):(2*pX+pXD)),size(XD0));
    VD = reshape(z((2*pX+pXD+1):end),size(XD0));
    for l=1:size(XD,1)
        XD(l,:) = max(1 - w(l)*lambda*mu./norm(XD(l,:)), 0).*XD(l,:);
        VD(l,:) = max(1 - w(l)*lambda*mu./norm(VD(l,:)), 0).*VD(l,:);
    end
    z_backward = [z(1:(2*pX));XD(:);VD(:)];
end
end