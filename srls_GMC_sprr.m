function [Xhat, Vhat, res_norm_hist] = srls_GMC_sprr(Y, X, lambda1, lambda2, varargin)

params = inputParser;
params.addParameter('gamma', 0.8, @(x) isnumeric(x));
params.addParameter('splitting', 'DY', @(x) ismember(x,{'DY'}));
params.addParameter('max_iter', 10000, @(x) isnumeric(x));
params.addParameter('tol_stop', 1e-5, @(x) isnumeric(x));
params.addParameter('early_termination', true, @(x) islogical(x));
params.addParameter('acceleration', 'aa2', @(x) ismember(x,{'original','aa2'}));
params.addParameter('mem_size', 5, @(x) isnumeric(x));
params.addParameter('eta', 1e-8, @(x) isnumeric(x));
params.parse(varargin{:});
gamma = params.Results.gamma;
splitting = params.Results.splitting;
max_iter = params.Results.max_iter;
tol_stop = params.Results.tol_stop;
early_termination = params.Results.early_termination;
acceleration = params.Results.acceleration;
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

n = size(X,1);
p = size(X,2);
L = size(Y,2);
Xt = X';
rho = norm(X)^2;
A = @(x) X*x;
AH = @(x) Xt*x;
XtY = Xt*Y;
mu = 1.99/(rho*(1-2*gamma+2*gamma^2)/(1-gamma));
zX0 = zeros(p,L);
zV0 = zeros(p,L);
z0 = [zX0(:);zV0(:)];
function z_proj = projection(z)
    zX = reshape(z(1:p*L),[p,L]);
    zV = reshape(z((p*L+1):end),[p,L]);
    zX = svt(zX, mu * lambda2);
    if gamma>0
        zV = svt(zV, mu * lambda2);
    end
    z_proj = [zX(:);zV(:)];
end
params_fixed.projection = @projection;
[z_lambda, iter, res_norm_hist] = fixed_iter(z0,@forward,@backward,params_fixed,acceleration);
Xhat = reshape(z_lambda(1:p*L),[p,L]);
Vhat = reshape(z_lambda((p*L+1):end),[p,L]);
fprintf('lambda1 = %f , lambda2 = %f solved in %d iterations\n', lambda1, lambda2, iter);

function z_forward = forward(z)
    X = reshape(z(1:p*L),[p,L]);
    V = reshape(z((p*L+1):end),[p,L]);
    zX = mu * ( AH(A(X + gamma*(V-X))) - XtY);
    if gamma>0
        zV = mu * ( gamma * AH(A(V-X)) );
    else
        zV = V;
    end
    z_forward = [zX(:);zV(:)];
end

function z_backward = backward(z)
    zX = reshape(z(1:p*L),[p,L]);
    zV = reshape(z((p*L+1):end),[p,L]);
    for i=1:p
        zX(i,:) = max(1 - lambda1*mu./norm(zX(i,:)), 0).*zX(i,:);
        if gamma>0
        	zV(i,:) = max(1 - lambda1*mu./norm(zV(i,:)), 0).*zV(i,:);
        end
    end
    z_backward = [zX(:);zV(:)];
end
end