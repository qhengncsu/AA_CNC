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
y = X*beta + sqrt(beta'*Sigma*beta/SNR)*randn(n,1);


% on a single lambda value, report the iterations
lambda_ratio = 0.1;
gamma = 0.8;

%% algorithm comparison on GMC: FB v.s. FB+AA
% vanilla
[x_fb, v_fb, res_norm_fb, t1] = srls_GMC_acc(y, X, lambda_ratio, 'type', 'single',...
                  'acceleration', 'original', "gamma", gamma, 'splitting', 'FB');


% AA
[x_fbaa, v_fbaa, res_norm_fbaa, t2] = srls_GMC_acc(y, X,  lambda_ratio, 'type','single',...
                  'acceleration', "aa2","gamma",gamma, 'splitting', 'FB');

%%
[x_fbf, v_fbf, res_norm_fbf, t3] = srls_GMC_acc(y, X, lambda_ratio, 'type', 'single',...
                  'acceleration', 'original',"gamma",gamma, 'splitting', 'FBF');


[x_fbfaa, v3_fbfaa, res_norm_fbfaa, t4] = srls_GMC_acc(y, X,  lambda_ratio, 'type','single',...
                  'acceleration', "aa2","gamma",gamma, 'splitting', 'FBF');

%%
%% algorithm comparison on GMC: DR v.s. DR+AA
% vanilla
[x_dr, v_dr, res_norm_dr, t5] = srls_GMC_acc(y, X, lambda_ratio, 'type', 'single',...
                  'acceleration', 'original',"gamma",gamma, 'splitting', 'DR');

% AA
[x_draa, v_draa, res_norm_draa, t6] = srls_GMC_acc(y, X,  lambda_ratio, 'type','single',...
                  'acceleration', "aa2","gamma",gamma, 'splitting', 'DR');

save("results/slr_single.mat","t1","t2","t3","t4","t5","t6","res_norm_fb","res_norm_fbaa",...
    "res_norm_fbf","res_norm_fbfaa","res_norm_dr","res_norm_draa");