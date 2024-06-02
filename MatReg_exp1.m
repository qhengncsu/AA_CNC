%% data generation
n = 500;
d1 = 64;
d2 = d1;
p = d1*d2;

rng(2024)
A = mvnrnd(zeros(p,1),eye(p), n);

X = zeros(d1);
row1 = 27;
row2 = 37;
X(row1:row2,10:54) = 1;
X(10:54,row1:row2) = 1;
imagesc(X);
colorbar;

% set parameters
gamma = 0.8;
lambda_seq = logspace(-1, 3, 20);

nReps = 20;
time_fb = zeros(nReps, 1);
time_fbaa = zeros(nReps, 1);
time_fbf = zeros(nReps, 1);
time_fbfaa = zeros(nReps, 1);

err_fb = zeros(nReps, 1);
err_fbaa = zeros(nReps, 1);
err_fbf = zeros(nReps, 1);
err_fbfaa = zeros(nReps, 1);

pred_fb = zeros(nReps, 1);
pred_fbaa = zeros(nReps, 1);
pred_fbf = zeros(nReps, 1);
pred_fbfaa = zeros(nReps, 1);

for k = 1:nReps
    y = A*X(:) + 0.5*randn(n,1);
    %% algorithm comparison on GMC: FB v.s. FB+AA
    % vanilla
    [X_fb, Xcell_fb, Vcell_fb, BICs_fb, t_fb] = GMC_matreg_path(y, A, d1, d2, 'lambda_seq',lambda_seq,...
                                    'acceleration', 'original', "gamma", gamma, 'splitting', 'FB');
   
    time_fb(k) = t_fb;
    err_fb(k) = norm(X_fb-X, 'fro')/norm(X, 'fro');
    pred_fb(k) = norm(y-A*X_fb(:), 'fro');
    
    % AA
    [X_fbaa, Xcell_fbaa, Vcell_fbaa, BICs_fbaa, t_fbaa] = GMC_matreg_path(y, A, d1, d2, 'lambda_seq',lambda_seq,...
                                          'acceleration', 'aa2', "gamma", gamma, 'splitting', 'FB');
    time_fbaa(k) = t_fbaa;
    err_fbaa(k) = norm(X_fbaa-X, 'fro')/norm(X, 'fro');
    pred_fbaa(k) = norm(y-A*X_fbaa(:), 'fro');
    
    %% algorithm comparison on GMC: FBF v.s. FBF+AA
    % vanilla
    [X_fbf, Xcell_fbf, Vcell_fbf, BICs_fbf, t_fbf] = GMC_matreg_path(y, A, d1, d2, 'lambda_seq',lambda_seq,...
                                         'acceleration', 'original', "gamma", gamma, 'splitting', 'FBF');
    
    time_fbf(k) = t_fbf;
    err_fbf(k) = norm(X_fbf-X, 'fro')/norm(X, 'fro');
    pred_fbf(k) = norm(y-A*X_fbf(:), 'fro');
    
    % AA
    [X_fbfaa, Xcell_fbfaa, Vcell_fbfaa, BICs_fbfaa, t_fbfaa] = GMC_matreg_path(y, A, d1, d2, 'lambda_seq',lambda_seq,...
                                          'acceleration', 'aa2', "gamma", gamma, 'splitting', 'FBF');
    time_fbfaa(k) = t_fbfaa;
    err_fbfaa(k) = norm(X_fbfaa-X, 'fro')/norm(X, 'fro');
    pred_fbfaa(k) = norm(y-A*X_fbfaa(:), 'fro');
    %%
    save("results/MatReg_exp1_n500.mat");
end

