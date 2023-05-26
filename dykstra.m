function [xhat, res_norm_hist] = dykstra(y, app, lambda, varargin)
% [xhat, vhat, res_norm_hist] = dykstra(y, app, lambda, varargin)
% INPUT
%   y 	    observed matrix
%   app     'closest kinship' or 'biclustering'
%   lambda  penalty parameter
% OUTPUT
%   xhat, vhat, res_norm_hist
params = inputParser;
params.addParameter('max_iter', 10000, @(x) isnumeric(x));
params.addParameter('tol_stop', 1e-5, @(x) isnumeric(x));
params.addParameter('early_termination', true, @(x) islogical(x));
params.addParameter('acceleration', 'aa2', @(x) ismember(x,{'original','aa2'}));
params.addParameter('mem_size', 5, @(x) isnumeric(x));
params.addParameter('eta', 1e-6, @(x) isnumeric(x));
params.addParameter('printevery', 100, @(x) isnumeric(x));
params.addParameter('z0', y(:), @(x) isnumeric(x));
params.parse(varargin{:});
max_iter = params.Results.max_iter;
tol_stop = params.Results.tol_stop;
early_termination = params.Results.early_termination;
acceleration = params.Results.acceleration;
mem_size = params.Results.mem_size;
eta = params.Results.eta;
printevery = params.Results.printevery;
z0 = params.Results.z0;
params_fixed = struct();
params_fixed.splitting = 'DK';
params_fixed.max_iter = max_iter;
params_fixed.tol = tol_stop;
params_fixed.early_termination = early_termination;
params_fixed.mem_size = mem_size;
params_fixed.verbose = true;
params_fixed.eta = eta;
params_fixed.printevery = printevery;
params_fixed.D = 1e6;
params_fixed.eta = 1e-6;
params_fixed.xi = 1e-14;
n1 = size(y,1);
n2 = size(y,2);
addpath(strcat(pwd,'/Solver.p/'));
addpath(strcat(pwd,'/utils/'));
[dim_row.d,dim_row.n] = size(y);
[dim_col.d,dim_col.n] = size(y');
k_n = 10;
phi = 0.01;
[row_weightVec,row_NodeArcMatrix] = compute_weight(y,k_n,phi,1);
[col_weightVec,col_NodeArcMatrix] = compute_weight(y',k_n,phi,1);
row_weightVec = row_weightVec/sum(row_weightVec,"all")/sqrt(n1);
col_weightVec = col_weightVec/sum(col_weightVec,"all")/sqrt(n2);

options.stoptol = 1e-6; %% tolerance for terminating the algorithm
options.num_k = k_n; %%number of nearest neighbors
A0_row = row_NodeArcMatrix;
Ainput_row.A = A0_row;
Ainput_row.Amap = @(x) x*A0_row;
Ainput_row.ATmap = @(x) x*A0_row';
Ainput_row.ATAmat = A0_row*A0_row'; %%graph Laplacian
Ainput_row.ATAmap = @(x) x*Ainput_row.ATAmat;
dim_row.E = length(row_weightVec);
A0_col = col_NodeArcMatrix;
Ainput_col.A = A0_col;
Ainput_col.Amap = @(x) x*A0_col;
Ainput_col.ATmap = @(x) x*A0_col';
Ainput_col.ATAmat = A0_col*A0_col'; %%graph Laplacian
Ainput_col.ATAmap = @(x) x*Ainput_col.ATAmat;
dim_col.E = length(col_weightVec);
options.use_kkt = 1;
options.printyes = 0;
options.printminoryes = 0;
options.maxiter = 1000;
options.admm_iter = 50;
[z_lambda, iter, res_norm_hist] = fixed_iter(-z0,@resolveP,@resolveQ,params_fixed,acceleration);
xhat = reshape(z_lambda(1:(n1*n2)),[n1,n2]);
fprintf('lambda = %f solved in %d iterations\n', lambda, iter);

function z_prox = resolveP(z)
    Z = reshape(z,[n1,n2]);
    if strcmp(app,'closest kinship')
        [V,D] = eig(Z);
        Z_prox = V*max(D,0)*V';
        z_prox = Z_prox(:);
    else
        evalc('[~,~,Z_prox,~,~,~] = SSNAL(Ainput_row,Z,dim_row,row_weightVec*lambda,options)');
        z_prox = Z_prox(:);
    end
end

function z_prox = resolveQ(z)
    Z = reshape(z,[n1,n2]);
    if strcmp(app,'closest kinship')
        Z_prox = max(Z,0);
        Z_prox = Z_prox - diag(diag(Z_prox)) + diag(0.5*ones(n1,1));
        z_prox = Z_prox(:);
    else
        Z = Z';
        evalc('[~,~,Z_prox,~,~,~] = SSNAL(Ainput_col,Z,dim_col,col_weightVec*lambda,options)');
        Z_prox = Z_prox';
        z_prox = Z_prox(:);
    end
end
end