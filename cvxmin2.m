function [xhat, res_norm_hist] = cvxmin2(Q, m, r, mu, z0, varargin)
% [xhat, res_norm_hist] = cvxmin2(Q, m, r, mu, z0, varargin)
% INPUT
%   Q 	    matrix in the objective function 
%   m       expected returns
%   r       expected total return 
%   mu      stepsize parameter
%   z0      initial guess
% OUTPUT
%   xhat
%   res_norm_hist
params = inputParser;
params.addParameter('max_iter', 10000, @(x) isnumeric(x));
params.addParameter('tol_stop', 1e-5, @(x) isnumeric(x));
params.addParameter('early_termination', true, @(x) islogical(x));
params.addParameter('splitting', 'DR', @(x) ismember(x,{'DR', 'DY'}));
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
params_fixed.projection = @projection;
p = length(z0);

if strcmp(splitting,'DR')
    % matrix for computing the proximal of quadratic term
    A_mat = sparse(Q + eye(p/3)/mu);
    % items for computing projection onto Ax=b
    AA = [eye(p/3) -eye(p/3) zeros(p/3);
           zeros(p/3) eye(p/3) -eye(p/3)];

    AA_pinv = pinv(AA);
    Afun = @(x) AA*x;
    Apfun = @(x) AA_pinv*x;
end



if strcmp(splitting,'DR')
    [xhat, iter, res_norm_hist] = fixed_iter(z0,@resolveP,@resolveQ,params_fixed,acceleration);
elseif strcmp(splitting,'DY')
    [xhat, iter, res_norm_hist] = fixed_iter(z0,@resolveP,@resolveQ,params_fixed,acceleration);
end

fprintf('Problem solved in %d iterations\n', iter);


% this is the forward
function z_prox = resolveP(z)
        
        if strcmp(splitting,'DR')
            z1_prox = projsplx(z(1:(p/3)));
            z2_prox = projhalf(z((p/3+1):(2*p/3)), -m, -r);
            
            b_vec = z((2*p/3+1):p)/mu;
            [z3_prox flag] = lsqr(A_mat, b_vec, 1e-10, 100);
            
            z_prox = [z1_prox; z2_prox; z3_prox];
        elseif strcmp(splitting,'DY')
            z_prox = z - mu*Q*z;
        end
end


% this is the backward
function z_prox = resolveQ(z)
    
    if strcmp(splitting,'DR')
        z_prox = z - Apfun(Afun(z));
    elseif strcmp(splitting,'DY')
        z_prox = projsplx(z);  %% simplex
    end
    
end


% this is the  projection
function z_proj = projection(z)
  
    z_proj = projhalf(z, -m ,-r);
    
end

end