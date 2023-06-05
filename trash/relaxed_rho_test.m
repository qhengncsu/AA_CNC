%% add path to the code 
addpath("Y:\MyDocuments\Xiaoqian\GMC-computation\code")

%% data generation
n = 500;
p = 50;
s = 5;
beta = [ones(s,1); zeros(p-s,1)];
rhoSeq = [0, 0.2, 0.5, 0.8];
numRho = length(rhoSeq);
SNR = 0.5;
nReps = 5;


time_GMC = nan(nReps, numRho);
time_Lasso = nan(nReps, numRho);

Err_rGMC = nan(nReps, numRho);
Err_GMC = nan(nReps, numRho);
Err_rLasso = nan(nReps, numRho);
Err_Lasso = nan(nReps, numRho);

Pd_rGMC = nan(nReps, numRho);
Pd_GMC = nan(nReps, numRho);
Pd_rLasso = nan(nReps, numRho);
Pd_Lasso = nan(nReps, numRho);

TP_rGMC = nan(nReps, numRho);
TP_GMC = nan(nReps, numRho);
TP_rLasso = nan(nReps, numRho);
TP_Lasso = nan(nReps, numRho);

FP_rGMC = nan(nReps, numRho);
FP_GMC = nan(nReps, numRho);
FP_rLasso = nan(nReps, numRho);
FP_Lasso = nan(nReps, numRho);

FN_rGMC = nan(nReps, numRho);
FN_GMC = nan(nReps, numRho);
FN_rLasso = nan(nReps, numRho);
FN_Lasso = nan(nReps, numRho);

F1_rGMC = nan(nReps, numRho);
F1_GMC = nan(nReps, numRho);
F1_rLasso = nan(nReps, numRho);
F1_Lasso = nan(nReps, numRho);

lambda_rGMC = nan(nReps, numRho);
lambda_GMC = nan(nReps, numRho);
lambda_rLasso = nan(nReps, numRho);
lambda_Lasso = nan(nReps, numRho);


Xlist = cell(numRho, 1);
ylist = cell(nReps, numRho);

%%
for j = 1: numRho
    rng(2023)
    % fix p
    rho = rhoSeq(j);
    Sigma = zeros(p, p);
    for k=1:p
        for l=1:p
            Sigma(k, l) = rho^(abs(k-l));
        end
    end

    X = mvnrnd(zeros(p,1),Sigma, n);
    
    Xlist{j, 1} = X;
    
    for i = 1:nReps
        
        rng(2023+i);
        y = X*beta + sqrt(beta'*Sigma*beta/SNR)*randn(n,1);
        ylist{i, j} = y;
     
   
       %% GMC (FB + AA2 + screning)
        alpha_seq = 0:0.25:1;
        t0 = tic;
        CV_GMC =  CV_relaxed_GMC(y, X, 'type', "single", 'gamma', 0.8, 'alpha_seq', alpha_seq,...
                          'lambda_min_ratio',0.005, 'acceleration', "aa2", "screen", true);
        t4 = toc(t0);      
        time_GMC(i, j) = t4;
        
        err_vec = zeros(length(alpha_seq),1);     
        for k=1:length(alpha_seq)
            err_vec(k, 1) = CV_GMC.min_err{k, 1};
        end
        [~, idx] = min(err_vec);
         best_beta = CV_GMC.bestbeta{idx, 1}';
         
        TP_rGMC(i, j) = length(intersect(find(CV_GMC.bestbeta{idx,1}~=0),find(beta~=0)));
        FP_rGMC(i, j) = length(intersect(find(CV_GMC.bestbeta{idx,1}~=0),find(beta==0)));
        FN_rGMC(i, j) = length(intersect(find(CV_GMC.bestbeta{idx,1}==0),find(beta~=0)));
        F1_rGMC(i, j) = 2*TP_rGMC(i, j) / (2*TP_rGMC(i, j) + FP_rGMC(i,j) + FN_rGMC(i, j));
        
        Err_rGMC(i, j) = norm(beta-best_beta, 'fro');
        Pd_rGMC(i, j) = norm(X*beta-X*best_beta, 'fro');
        lambda_rGMC(i, j) = CV_GMC.min_lambda{idx, 1};
        
        
        % GMC
        TP_GMC(i, j) = length(intersect(find(CV_GMC.bestbeta{end,1}~=0),find(beta~=0)));
        FP_GMC(i, j) = length(intersect(find(CV_GMC.bestbeta{end,1}~=0),find(beta==0)));
        FN_GMC(i, j) = length(intersect(find(CV_GMC.bestbeta{end,1}==0),find(beta~=0)));
        F1_GMC(i, j) = 2*TP_GMC(i, j) / (2*TP_GMC(i, j) + FP_GMC(i,j) + FN_GMC(i, j));
        
        Err_GMC(i, j) = norm(beta-CV_GMC.bestbeta{end,1}', 'fro');
        Pd_GMC(i, j) = norm(X*beta-X*CV_GMC.bestbeta{end,1}', 'fro');
        lambda_rGMC(i, j) = CV_GMC.min_lambda{end, 1};
   
        
        %% GMC (FB + AA2 + screning)
        alpha_seq = 0:0.25:1;
        t0 = tic;
        CV_Lasso =  CV_relaxed_GMC(y, X, 'type', "single", 'gamma', 0, 'alpha_seq', alpha_seq,...
                          'lambda_min_ratio',0.005, 'acceleration', "aa2", "screen", true);
        t4 = toc(t0);      
        time_Lasso(i, j) = t4;
        
        err_vec = zeros(length(alpha_seq),1);     
        for k=1:length(alpha_seq)
            err_vec(k, 1) = CV_Lasso.min_err{k, 1};
        end
        [~, idx] = min(err_vec);
         best_beta = CV_Lasso.bestbeta{idx, 1}';
         
        TP_rLasso(i, j) = length(intersect(find(CV_Lasso.bestbeta{idx,1}~=0),find(beta~=0)));
        FP_rLasso(i, j) = length(intersect(find(CV_Lasso.bestbeta{idx,1}~=0),find(beta==0)));
        FN_rLasso(i, j) = length(intersect(find(CV_Lasso.bestbeta{idx,1}==0),find(beta~=0)));
        F1_rLasso(i, j) = 2*TP_rLasso(i, j) / (2*TP_rLasso(i, j) + FP_rLasso(i,j) + FN_rLasso(i, j));
        
        Err_rLasso(i, j) = norm(beta-best_beta, 'fro');
        Pd_rLasso(i, j) = norm(X*beta-X*best_beta, 'fro');
        lambda_rLasso(i, j) = CV_Lasso.min_lambda{idx, 1};
        
        
        % Lasso
        TP_Lasso(i, j) = length(intersect(find(CV_Lasso.bestbeta{end,1}~=0),find(beta~=0)));
        FP_Lasso(i, j) = length(intersect(find(CV_Lasso.bestbeta{end,1}~=0),find(beta==0)));
        FN_Lasso(i, j) = length(intersect(find(CV_Lasso.bestbeta{end,1}==0),find(beta~=0)));
        F1_Lasso(i, j) = 2*TP_Lasso(i, j) / (2*TP_Lasso(i, j) + FP_Lasso(i,j) + FN_Lasso(i, j));
        
        Err_Lasso(i, j) = norm(beta-CV_Lasso.bestbeta{end,1}', 'fro');
        Pd_Lasso(i, j) = norm(X*beta-X*CV_Lasso.bestbeta{end,1}', 'fro');
        lambda_rLasso(i, j) = CV_Lasso.min_lambda{end, 1};
       
        %% 
       save relaxed_rho_test.mat
    end
    
       save relaxed_rho_test.mat
   
    
end