function [xhat, res_norm_hist] = cvx_min(A, b, mu, z0, app, varargin)
% [xhat, res_norm_hist] = cvxmin(A, b, varargin)
% INPUT
%   A 	    matrix in the objective function
%   b       vector in the objective function; Y for SRRR
%   mu      stepsize parameter
%   z0      initial guess; vectoried X for SRRR
%   app     two options: 'NNLS' or 'SRRR'
% OUTPUT
%   xhat
%   res_norm_hist
params = inputParser;
params.addParameter('lambda1', 0, @(x) isnumeric(x));
params.addParameter('lambda2', 0, @(x) isnumeric(x));
params.addParameter('max_iter', 10000, @(x) isnumeric(x));
params.addParameter('tol_stop', 1e-5, @(x) isnumeric(x));
params.addParameter('early_termination', true, @(x) islogical(x));
params.addParameter('splitting', 'DR', @(x) ismember(x,{'DR','FB', 'DY'}));
params.addParameter('acceleration', 'aa2', @(x) ismember(x,{'original','aa2'}));
params.addParameter('mem_size', 10, @(x) isnumeric(x));
params.addParameter('eta', 1e-8, @(x) isnumeric(x));
params.addParameter('printevery', 100, @(x) isnumeric(x));
params.parse(varargin{:});

% load values of those argumens
max_iter = params.Results.max_iter;
tol_stop = params.Results.tol_stop;
early_termination = params.Results.early_termination;
splitting = params.Results.splitting;
acceleration = params.Results.acceleration;
mem_size = params.Results.mem_size;
eta = params.Results.eta;
printevery = params.Results.printevery;
lambda1 = params.Results.lambda1;
lambda2 = params.Results.lambda2;

% set parameters for fixed iteration
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
params_fixed.D = 10;
params_fixed.xi = 1e-14;
params_fixed.projection = @projection;

% compute required iterms for different applications
if strcmp(app,'NNLS')
    p = size(A,2);
    if strcmp(splitting,'DR')
        A_mat = sparse(A'*A + eye(p)/mu);
    elseif strcmp(splitting,'FB')
        Afun = @(x) A*x;
        AH = @(x) A'*x;
    end
    Atb = A'*b;
elseif strcmp(app,'SRRR')
    [n, p] = size(A);
    L = length(b)/n;
    Y = reshape(b, [n, L]);
    
    rho = norm(A)^2;
    mu = 1.99/rho;
    
    if strcmp(splitting,'DR')
        A_mat = sparse(A'*A + eye(p)/mu);
        A_inv = inv(A_mat);
        % items for computing projection onto Ax=b
        AA = sparse([eye(p*L) -eye(p*L) zeros(p*L);
            zeros(p*L) eye(p*L) -eye(p*L)]);
        AAT = AA*AA';
    elseif strcmp(splitting,'DY')
        Afun = @(x) A*x;
        AH = @(x) A'*x;
    end
    AtY = A'*Y;
end


[xhat, iter, res_norm_hist] = fixed_iter(z0,@resolveP,@resolveQ,params_fixed,acceleration);

fprintf('Problem solved in %d iterations\n', iter);


% this is the forward step or resolvent of P
function z_prox = resolveP(z)
    if strcmp(app,'NNLS')
        if strcmp(splitting,'DR')
            b_vec = Atb + z/mu;
            [z_prox, ~]= lsqr(A_mat, b_vec);
        elseif strcmp(splitting,'FB')
            z_prox = z - mu * (AH(Afun(z)) - Atb);
        end
    elseif strcmp(app,'SRRR')       
        if strcmp(splitting,'DR')
            % convert the vector to matrix
            X1 = reshape(z(1:(p*L)), [p, L]);
            X2 = reshape(z((p*L+1):(2*p*L)), [p, L]);
            X3 = reshape(z((2*p*L+1):(3*p*L)), [p, L]);
            % update X1
            zX1 =  A_inv*(AtY + X1/mu);
            zX1_prox = zX1(:);
            % update X2
            zX2 = proj_L2(X2, mu * lambda1);
            zX2_prox = zX2(:);
            % update X3
            zX3 = svt(X3, mu * lambda2);
            zX3_prox = zX3(:);
            % combine all updates
            z_prox = [zX1_prox; zX2_prox; zX3_prox];
        elseif strcmp(splitting,'DY')    
            X = reshape(z,[p,L]);
            zX = X- mu * ( AH(Afun(X)) - AtY);
            z_prox = zX(:);
        end
    end
end

% this is backward step or resolvent of Q
function z_prox = resolveQ(z)
    if strcmp(app,'NNLS')
        z_prox = max(z, 0);
    elseif strcmp(app,'SRRR')
        if strcmp(splitting,'DR')
            z_prox = z - AA'*(AAT\(AA*z));
        elseif strcmp(splitting,'DY')
            X = reshape(z,[p,L]);
            zX = proj_L2(X, mu * lambda1);
            z_prox = zX(:);
        end
    end
end

% this is projection step in DY
function z_proj = projection(z)
    X = reshape(z,[p,L]);
    zX = svt(X, mu * lambda2);
    z_proj = zX(:);
end

% compute the proximal of row-wise sparsity penalty
function X_proj = proj_L2(X, t)
    r = size(X, 1);
    X_proj = X;
    for i=1:r
        X_proj(i,:) = max(1 - t./norm(X(i,:)), 0).*X(i,:);
    end
end

end