%% add path to the code 
%addpath("Y:\MyDocuments\Xiaoqian\GMC-computation\code")
%% data generation
n = 2000;
p = 10000;
SNR = 1;
rng(2023)
rho = 0.3;
Sigma = zeros(p, p);
for k=1:p
    for l=1:p
        Sigma(k, l) = rho^(abs(k-l));
    end
end
X = mvnrnd(zeros(p,1),Sigma, n);

gp_size = 50;
gp_num = p/gp_size;

beta = [ones(gp_size,1); -ones(gp_size,1);  zeros(gp_size*(gp_num-2),1)];

groups = cell(gp_num,1);
for i=1:gp_num
    groups{i} = ((i-1)*gp_size+1):(i*gp_size);
end

%% Outputs
nReps = 20;
ylist = cell(nReps, 1);

TP_GMC = nan(nReps, 1);
FP_GMC = nan(nReps, 1);
FN_GMC = nan(nReps, 1);
F1_GMC = nan(nReps, 1);
Err_GMC = nan(nReps, 1);
Pd_GMC = nan(nReps, 1);
lambda_GMC = nan(nReps, 1);

TP_grGMC = nan(nReps, 1);
FP_grGMC = nan(nReps, 1);
FN_grGMC = nan(nReps, 1);
F1_grGMC = nan(nReps, 1);
Err_grGMC = nan(nReps, 1);
Pd_grGMC = nan(nReps, 1);
lambda_grGMC = nan(nReps, 1);

%% For loop
for i = 1: nReps
    
    rng(2023+i);
    y = X*beta + sqrt(beta'*Sigma*beta/SNR)*randn(n,1);
    ylist{i, 1} = y;
    
    %% GMC
    t0 = tic;
    cv_GMC = CV_GMC_path(y, X, 'type', "single", 'gamma', 0.8,...
                          'lambda_min_ratio',0.01, 'acceleration', "aa2", "screen", false);
    t1 = toc(t0);
    
    
    best_beta_GMC = cv_GMC.bestbeta';
    
    TP_GMC(i, 1) = length(intersect(find(best_beta_GMC~=0),find(beta~=0)));
    FP_GMC(i, 1) = length(intersect(find(best_beta_GMC~=0),find(beta==0)));
    FN_GMC(i, 1) = length(intersect(find(best_beta_GMC==0),find(beta~=0)));
    F1_GMC(i, 1) = 2*TP_GMC(i, 1) / (2*TP_GMC(i, 1) + FP_GMC(i, 1) + FN_GMC(i, 1));
    
    Err_GMC(i, 1) = norm(beta-best_beta_GMC, 'fro');
    Pd_GMC(i, 1) = norm(X*beta-X*best_beta_GMC, 'fro');
    lambda_GMC(i, 1) = cv_GMC.min_lambda;
    
    
    save GMC_cv.mat
    %% grGMC
    t0 = tic;
    cv_grGMC = CV_GMC_path(y, X, 'type', "group", 'groups', groups, 'gamma', 0.8,...
                            'lambda_min_ratio',0.01, 'acceleration', "aa2", "screen", false);
    t3 = toc(t0);
    
    
    best_beta_grGMC = cv_grGMC.bestbeta';
    
    TP_grGMC(i, 1) = length(intersect(find(best_beta_grGMC~=0),find(beta~=0)));
    FP_grGMC(i, 1) = length(intersect(find(best_beta_grGMC~=0),find(beta==0)));
    FN_grGMC(i, 1) = length(intersect(find(best_beta_grGMC==0),find(beta~=0)));
    F1_grGMC(i, 1) = 2*TP_grGMC(i, 1) / (2*TP_grGMC(i, 1) + FP_grGMC(i, 1) + FN_grGMC(i, 1));
    
    Err_grGMC(i, 1) = norm(beta-best_beta_grGMC, 'fro');
    Pd_grGMC(i, 1) = norm(X*beta-X*best_beta_grGMC, 'fro');
    lambda_grGMC(i, 1) = cv_grGMC.min_lambda;
    
    save GMC_cv.mat
    
    
end
