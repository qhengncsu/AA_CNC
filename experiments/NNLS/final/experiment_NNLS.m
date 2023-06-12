
%%
n = 3e4;
p = 3e4;

%%
rng(2023)
% use sparse matrix
s = 0.5; % sparsity ratio, a alrger sparsity ratio leads to higher speed up ratio
A = randn(n, p);
a = A(:);
a(randsample(n*p, s*n*p)) = 0;
A = sparse(reshape(a, [n, p]));

b = randn(n, 1);
tic;
mu = 1.99/normest(A)^2;
toc;
%mu = 0.1;
%z0 = lsqr(A, b, 1e-10, 100);
z0 = zeros(p, 1);
app = 'NNLS';

%% DR
t0 = tic;
[xhat_DR_aa, res_norm_hist_DR_aa] = cvxmin(A, b, mu, z0, app, 'max_iter', 1e4);
t_DR_aa = toc(t0);

t0 = tic;
[xhat_DR, res_norm_hist_DR] = cvxmin(A, b, mu, z0, app, 'max_iter', 1e4, 'acceleration', 'original');
t_DR = toc(t0);

%% evaluation of DR
obj_DR_aa = norm(A*xhat_DR_aa - b, 'fro')

obj_DR = norm(A*xhat_DR - b, 'fro')

err_DR_aa = norm(min(xhat_DR_aa, 0))^2
err_DR = norm(min(xhat_DR, 0))^2


%% FB
t0 = tic;
[xhat_FB_aa, res_norm_hist_FB_aa] = cvxmin(A, b, mu, z0, app, 'max_iter', 1e4, 'splitting', 'FB');
t_FB_aa = toc(t0);

t0 = tic;
[xhat_FB, res_norm_hist_FB] = cvxmin(A, b, mu, z0, app, 'max_iter', 1e4, 'acceleration', 'original', 'splitting', 'FB');
t_FB = toc(t0);

%% evaluation of FB
obj_FB_aa = norm(A*xhat_FB_aa - b, 'fro')

obj_FB = norm(A*xhat_FB - b, 'fro')

err_FB_aa = norm(min(xhat_FB_aa, 0))^2
err_FB = norm(min(xhat_FB, 0))^2


%%
 save NNLS_results_n3e4_s5.mat
