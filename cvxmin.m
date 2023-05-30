function [xhat, res_norm_hist] = cvxmin(A, b, mu, z0, app, varargin)
% [xhat, res_norm_hist] = cvxmin(A, b, varargin)
% INPUT
%   A 	    matrix in the objective function 
%   b       vector in the objective function
%   app     'NNLS' or ''
%   mu      stepsize parameter
% OUTPUT
%   xhat
%   res_norm_hist
params = inputParser;
params.addParameter('max_iter', 10000, @(x) isnumeric(x));
params.addParameter('tol_stop', 1e-5, @(x) isnumeric(x));
params.addParameter('early_termination', true, @(x) islogical(x));
params.addParameter('splitting', 'DR', @(x) ismember(x,{'DR','FB', 'FBF', 'DY', 'DK'}));
params.addParameter('acceleration', 'aa2', @(x) ismember(x,{'original','aa2'}));
params.addParameter('mem_size', 10, @(x) isnumeric(x));
params.addParameter('eta', 1e-8, @(x) isnumeric(x));
params.addParameter('printevery', 100, @(x) isnumeric(x));
params.parse(varargin{:});

max_iter = params.Results.max_iter;
tol_stop = params.Results.tol_stop;
early_termination = params.Results.early_termination;
splitting = params.Results.splitting;
acceleration = params.Results.acceleration;
mem_size = params.Results.mem_size;
eta = params.Results.eta;
printevery = params.Results.printevery;

params_fixed = struct();
params_fixed.max_iter = max_iter;
params_fixed.tol = tol_stop;
params_fixed.early_termination = early_termination;
params_fixed.splitting = splitting;
params_fixed.acceleration = acceleration;
params_fixed.mem_size = mem_size;
params_fixed.verbose = true;
params_fixed.eta = eta;
params_fixed.printevery = printevery;
params_fixed.D = 1;
params_fixed.eta = eta;
params_fixed.xi = 1e-14;
p = size(A,2);
if strcmp(splitting,'DR')
    A_mat = sparse(A'*A + eye(p)/mu);
else
    Afun = @(x) A*x;
    AH = @(x) A'*x;
end
Atb = A'*b;
[xhat, iter, res_norm_hist] = fixed_iter(z0,@resolveP,@resolveQ,params_fixed,acceleration);

fprintf('Problem solved in %d iterations\n', iter);


% this is the proximal of mu*f
function z_prox = resolveP(z)
    if strcmp(app,'NNLS')
        if strcmp(splitting,'DR')
            b_vec = Atb + z/mu;
            [z_prox, flag]= lsqr(A_mat, b_vec, 1e-10, 100);
        elseif strcmp(splitting,'FB')
            z_prox = z - mu * (AH(Afun(z)) - Atb);
        end
    else
    end
end


% this is the proximal of mu*g
function z_prox = resolveQ(z)
    
    if strcmp(app,'NNLS')
        z_prox = max(z, 0);
    else
    end
    
end


end
