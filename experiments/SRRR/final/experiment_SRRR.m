
%%
load data_mice
%%
A = normalize(X);
Y = normalize(Y);
app = 'SRRR';
%%
p = size(A, 2); % A is n-by-p
L = size(Y, 2); % Y is n-by-L
b = Y(:);
mu = 0; % give a random mu for SPRR, the algorithm will compute a good one
lambda1 = 46;
lambda2 = 129;
%% DY
z0 = zeros(p*L, 1);

t0 = tic;
[xhat_DY, res_norm_hist_DY] = cvxmin(A, b, mu, z0, app, 'lambda1', lambda1, 'lambda2', lambda2,...
                                 'max_iter', 1e4, 'acceleration', 'original', 'splitting', 'DY');
t_DY = toc(t0);


t0 = tic;
[xhat_DY_aa, res_norm_hist_DY_aa] = cvxmin(A, b, mu, z0, app, 'lambda1', lambda1, 'lambda2', lambda2,...
                                     'max_iter', 1e4, 'splitting', 'DY');
t_DY_aa = toc(t0);

%% see of objective values from DY
X_DY = reshape(xhat_DY, [p, L]);
obj_DY = obj_SPRR(Y, A, X_DY, lambda1, lambda2)

X_DY_aa = reshape(xhat_DY_aa, [p, L]);
obj_DY_aa = obj_SPRR(Y, A, X_DY_aa, lambda1, lambda2)

%%

%% DR
z0 = zeros(3*p*L, 1);
mu = 1;
t0 = tic;
[xhat_DR, res_norm_hist_DR] = cvxmin(A, b, mu, z0, app, 'lambda1', lambda1, 'lambda2', lambda2,...
                                 'max_iter', 1e4, 'acceleration', 'original', 'splitting', 'DR');
t_DR = toc(t0);


t0 = tic;
[xhat_DR_aa, res_norm_hist_DR_aa] = cvxmin(A, b, mu, z0, app, 'lambda1', lambda1, 'lambda2', lambda2,...
                                     'max_iter', 1e4, 'splitting', 'DR');
t_DR_aa = toc(t0);

%% see of objective values from DR
z = ( xhat_DR(1:(p*L)) + xhat_DR((p*L+1):(2*p*L)) + xhat_DR((2*p*L+1):end))/3;
X_DR = reshape(z, [p, L]);
obj_DR = obj_SPRR(Y, A, X_DR, lambda1, lambda2)

z_aa = ( xhat_DR_aa(1:(p*L)) + xhat_DR_aa((p*L+1):(2*p*L)) + xhat_DR_aa((2*p*L+1):end))/3;
X_DR_aa = reshape(z_aa, [p, L]);
obj_DR_aa = obj_SPRR(Y, A, X_DR_aa, lambda1, lambda2)


%%

save SRRR.mat
