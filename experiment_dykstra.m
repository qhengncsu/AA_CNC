data_dir = pwd;
addpath(strcat(pwd,'/Solver.p/'));
addpath(strcat(pwd,'/utils/'));
dataMatrix = csvread('lung500.csv',1,1);
%dataMatrix = checkerboard(50,2,2)+0.5*randn(200,200)*std(checkerboard(50,2,2),[],'all');
global dim_row dim_col d n gamma_row gamma_col row_weightVec1 col_weightVec1 Ainput_row Ainput_col options;
[dim_row.d,dim_row.n] = size(dataMatrix);
[dim_col.d,dim_col.n] = size(dataMatrix');
d = dim_row.d;
n = dim_row.n;
k_n = 10;
phi = 0.01;
[row_weightVec1,row_NodeArcMatrix] = compute_weight(dataMatrix,k_n,phi,1);
[col_weightVec1,col_NodeArcMatrix] = compute_weight(dataMatrix',k_n,phi,1);
row_weightVec1 = row_weightVec1/sum(row_weightVec1,"all")*1/sqrt(d);
col_weightVec1 = col_weightVec1/sum(col_weightVec1,"all")*1/sqrt(n);

options.stoptol = 1e-6; %% tolerance for terminating the algorithm
options.num_k = k_n; %%number of nearest neighbors
A0_row = row_NodeArcMatrix;
Ainput_row.A = A0_row;
Ainput_row.Amap = @(x) x*A0_row;
Ainput_row.ATmap = @(x) x*A0_row';
Ainput_row.ATAmat = A0_row*A0_row'; %%graph Laplacian
Ainput_row.ATAmap = @(x) x*Ainput_row.ATAmat;
dim_row.E = length(row_weightVec1);
A0_col = col_NodeArcMatrix;
Ainput_col.A = A0_col;
Ainput_col.Amap = @(x) x*A0_col;
Ainput_col.ATmap = @(x) x*A0_col';
Ainput_col.ATAmat = A0_col*A0_col'; %%graph Laplacian
Ainput_col.ATAmap = @(x) x*Ainput_col.ATAmat;
dim_col.E = length(col_weightVec1);
options.use_kkt = 1;
options.printyes = 0;
options.printminoryes = 0;
options.maxiter = 1000;
options.admm_iter = 50;
gamma_row = 3e5;
gamma_col = 3e5;
params_fixed.splitting = 'DK';
params_fixed.max_iter = 100;
params_fixed.tol = 1e-5;
params_fixed.early_termination = true;
params_fixed.mem_size = 5;
params_fixed.verbose = true;
params_fixed.eta = 1e-8;
params_fixed.printevery = 10;
z0 = -dataMatrix(:);
[z_lambda, iter, res_norm_hist] = fixed_iter(z0,@forward,@backward,params_fixed,'aa2');
xhat1 = reshape(z_lambda,[d,n]);
imshow(xhat1')
tolClustering = 1e4*options.stoptol;
[cluster_id, num_cluster] = find_cluster(xhat1,tolClustering);

function z_prox = forward(z)
    global dim_row d n gamma_row row_weightVec1 Ainput_row options;
    X = reshape(z,[d,n]);
    weightVec = gamma_row*row_weightVec1;
    evalc('[~,~,X_prox,~,~,~] = SSNAL(Ainput_row,X,dim_row,weightVec,options)');
    z_prox = X_prox(:);
end

function z_prox = backward(z)
    global dim_col d n gamma_col col_weightVec1 Ainput_col options;
    X = reshape(z,[d,n]);
    X = X';
    weightVec = gamma_col*col_weightVec1;
    evalc('[~,~,X_prox,~,~,~] = SSNAL(Ainput_col,X,dim_col,weightVec,options)');
    X_prox = X_prox';
    z_prox = X_prox(:);
end

