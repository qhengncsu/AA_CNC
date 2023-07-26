% %% add path to the code 
% addpath('/rsrch6/home/bcb/xliu31/GMC/code')
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

%% Output
nReps = 20;
TM_sg_original = nan(nReps, 1);
TM_sg_aa2 = nan(nReps, 1);
diff_sg = nan(nReps, 1);


TM_gp_original = nan(nReps, 1);
TM_gp_aa2 = nan(nReps, 1);
diff_gp = nan(nReps, 1);

%%
for i = 1: nReps
   
    
    rng(2023+i);
    y = X*beta + sqrt(beta'*Sigma*beta/SNR)*randn(n,1);
    
    
    %% GMC (FB original)
    t0 = tic;
    [xhat1_matrix, vhat1_matrix, intercept1, lambda_seq1] = srls_GMC_path(y, X, 'type', "single",...
                                        'gamma', 0.8,'lambda_min_ratio',0.01, 'acceleration', "original", "screen",false);
    t1 = toc(t0);
    TM_sg_original(i, 1) = t1;
 
    
    % GMC (FB aa2)
    t0 = tic;
    [xhat2_matrix, vhat2_matrix, intercept2, lambda_seq2] = srls_GMC_path(y, X, 'type', "single",...
                                        'gamma', 0.8,'lambda_min_ratio',0.01, 'acceleration', "aa2", "screen",false);
    t2 = toc(t0);
    TM_sg_aa2(i, 1) = t2;
    
    %
    diff_sg(i, 1) = norm(xhat1_matrix - xhat2_matrix, 'fro');
    
    %% grGMC (FB original)
    t0 = tic;
    [xhat3_matrix, vhat3_matrix, intercept3, lambda_seq3] = srls_GMC_path(y, X, 'type', "group", 'groups',groups,...
                                    'gamma', 0.8,'lambda_min_ratio',0.01, 'acceleration', "original", "screen",false);
    t3 = toc(t0);
    TM_gp_original(i, 1) = t3;
 
    
    % grGMC (FB aa2)
    t0 = tic;
    [xhat4_matrix, vhat4_matrix, intercept4, lambda_seq4] = srls_GMC_path(y, X, 'type', "group", 'groups',groups,...
                                     'gamma', 0.8,'lambda_min_ratio',0.01, 'acceleration', "aa2", "screen",false);
    t4 = toc(t0);
    TM_gp_aa2(i, 1) = t4;
 
    %
    diff_gp(i, 1) = norm(xhat3_matrix - xhat4_matrix);
    
    
end
save("results/GMC_path_time.mat","TM_sg_original","TM_sg_aa2","TM_gp_original","TM_gp_aa2")

