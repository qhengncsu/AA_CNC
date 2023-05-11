addpath(strcat(pwd,'/Solver.p/'));
addpath(strcat(pwd,'/utils/'));
global mu AHy n p row_edges col_edges AT A AH lambda_row lambda_col row_weightVec col_weightVec;
%Y = checkerboard(50,2,2)+randn(200,200)*std(checkerboard(50,2,2),[],'all');
Y = csvread('lung500.csv',1,1);
Y = Y';
[n,p] = size(Y);
k_n = 10;
phi = 0.01;
[row_weightVec,row_NodeArcMatrix] = compute_weight(Y',k_n,phi,1);
[col_weightVec,col_NodeArcMatrix] = compute_weight(Y,k_n,phi,1);
row_weightVec = row_weightVec/sum(row_weightVec,"all")/sqrt(p);
col_weightVec = col_weightVec/sum(col_weightVec,"all")/sqrt(n);
row_NodeArcMatrix = row_NodeArcMatrix';
col_NodeArcMatrix = col_NodeArcMatrix';
Ar = kron(speye(p),row_NodeArcMatrix);
Ac = kron(col_NodeArcMatrix,speye(n));
row_edges = length(row_weightVec);
col_edges = length(col_weightVec);
AT = [Ar;Ac]';
A = @(x) AT*x;
AH = @(x) AT'*x;
AHy = AH(Y(:));
rho = normest(AT)^2;
mu = 1/rho;
lambda_row = 3e5;
lambda_col = 3e5;
params_fixed.splitting = 'FB';
params_fixed.max_iter = 20000;
params_fixed.tol = 1e-5;
params_fixed.early_termination = true;
params_fixed.mem_size = 2;
params_fixed.verbose = true;
params_fixed.eta = 1e-4;
params_fixed.printevery = 100;
z0 = zeros(size(AT,2),1);
% tic
% [z_lambda1, iter1, res_norm_hist1] = fixed_iter(z0,@forward,@backward,params_fixed,'nesterov');
% toc
tic
[z_lambda2, iter2, res_norm_hist2] = fixed_iter(z0,@forward,@backward,params_fixed,'aa2');
toc
X = reshape(Y(:) - A(z_lambda2),[n,p]);
imshow(X)
tolClustering = 1e-2;
[cluster_id, num_cluster] = find_cluster(X',tolClustering);

opts = [];
opts.tol = 1e-4;
opts.recordObjective = true; %  Record the objective function so we can plot it
opts.verbose = true;
opts.stringHeader='    ';
b = Y(:);
f    = @(z) .5*norm(b-z,'fro')^2;
grad = @(z) z-b;
g = @(x) 0;
% proxg(z,t) = argmin t*mu*|x|+.5||x-z||^2
prox = @(z,t) backward(z);

%% Call solver
tic
[solution, outs] = fasta(A,AH,f,grad,g,prox,z0,opts);
toc


function z_forward = forward(z)
    global mu AHy A AH;
    z_forward = z - mu * ( AH(A(z)) - AHy);
end

function z_backward = backward(z)
    global mu AHy n p row_edges col_edges AT A AH lambda_row lambda_col row_weightVec col_weightVec;
    Z_row = reshape(z(1:(row_edges*p)),[row_edges,p]);
    Z_col = reshape(z((row_edges*p+1):size(AT,2)),[n,col_edges]);
    Z_row_prox = proj_l2(Z_row', lambda_row*row_weightVec)';
    Z_col_prox = proj_l2(Z_col, lambda_col*col_weightVec);
    % row_norms = vecnorm(Z_row);
    % col_norms = vecnorm(Z_col);
    % Z_row_scale = min(lambda_row*row_weightVec./row_norms, 1);
    % Z_row_prox = Z_row.*Z_row_scale;
    % Z_col_scale = min(lambda_col*col_weightVec./col_norms, 1);
    % Z_col_prox = Z_col.*Z_col_scale;
    z_backward = [Z_row_prox(:);Z_col_prox(:)];
end

