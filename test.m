%% add path to the code 
addpath("Y:\MyDocuments\Xiaoqian\GMC-computation\code")

%% data generation
n = 1000;
p = 5000;
SNR = 1;

X = 5*randn(n,p);
beta = [ones(p/100,1); -2*ones(p/100, 1); zeros(p*49/50,1)];
y = X*beta + randn(n,1)*std(X*beta)/SNR;

groups = cell(p/100,1);
for i=1:p/100
    groups{i} = ((i-1)*10+1):(i*10);
end
% data standarization
% center = mean(X);
% scale = sqrt(sum((X - center).^2)/n);
% XX = (X - center)./scale;
% yy = y -mean(y);
% yy = (y-mean(y))/std(y);
% XX = normc(X);  if normalize this way, the nonzero elements are different
% with/without normalization
%% algorithm comparison: vanilla, Nesterov, AA on GMC
% on a single lambda value, report the iterations
lambda_ratio = 0.1 ;
gamma = 0.8;
% vanilla
t0 = tic;
[x1_sg, v1_sg, res_norm_sg1] = srls_GMC_acc(y, X, lambda_ratio, 'type', 'single', 'acceleration', 'original',"gamma",gamma);
t1 = toc(t0);
% x1 = x1_sg./scale';
% Nesterov
t0 = tic;
[x2_sg, v2_sg, res_norm_sg2] = srls_GMC_acc(y, X, lambda_ratio, 'type', "single", 'acceleration', "inertia","gamma",gamma);
t2 = toc(t0);
% AA
t0 = tic;
[x3_sg, v3_sg, res_norm_sg3] = srls_GMC_acc(y, X,  lambda_ratio, 'type', "single", 'acceleration', "aa2","gamma",gamma);
t3 = toc(t0);

%% Plots
figure;
plot(log(res_norm_sg1), 'k-', 'LineWidth', 1.5)
hold on; 
plot(log(res_norm_sg2), 'b-', 'LineWidth', 1.5)
hold on;
plot(log(res_norm_sg3), 'r-',  'LineWidth',1.5)
%xlabel('$iteration$','Interpreter','latex');
xlabel('Iteration');
ylabel('Norm of redidual');
legend(['Vanilla FB (t=', num2str(round(t1, 2)), ')'],...
       ['Nesterov+FB (t=', num2str(round(t2, 2)), ')'],... 
       ['AA2+FB (t=', num2str(round(t3, 2)), ')'], ...
       'Location', 'best');
       %'Location', 'best','Interpreter','latex');
title(['Single solution of GMC at lambda\_ratio = ', num2str(lambda_ratio)])
ax = gca; 
ax.FontSize = 12; 



%% algorithm comparison: vanilla, Nesterov, AA on group GMC
% on a single lambda value, report the iterations
lambda_ratio = 0.1;
gamma = 0.8;
% vanilla
t0 = tic;
[x1_gp, v1_gp, res_norm_gp1] = srls_GMC_acc(y, X, lambda_ratio, 'type', 'grouped', 'groups', groups, 'acceleration', 'original',"gamma",gamma);
t1 = toc(t0);
% Nesterov
t0 = tic;
[x2_gp, v2_gp, res_norm_gp2] = srls_GMC_acc(y, X, lambda_ratio, 'type', "grouped", 'groups', groups,  'acceleration', "inertia","gamma",gamma);
t2 = toc(t0);
% AA
t0 = tic;
[x3_gp, v3_gp, res_norm_gp3] = srls_GMC_acc(y, X,  lambda_ratio, 'type', "grouped", 'groups', groups, 'acceleration', "aa2","gamma",gamma);
t3 = toc(t0);

%% Plots
figure;
plot(log(res_norm_gp1), 'k-', 'LineWidth', 1.5)
hold on; 
plot(log(res_norm_gp2), 'b-', 'LineWidth', 1.5)
hold on;
plot(log(res_norm_gp3), 'r-',  'LineWidth',1.5)
%xlabel('$iteration$','Interpreter','latex');
xlabel('Iteration');
ylabel('Norm of redidual');
legend(['Vanilla FB (t=', num2str(round(t1, 2)), ')'],...
       ['Nesterov+FB (t=', num2str(round(t2, 2)), ')'],... 
       ['AA2+FB (t=', num2str(round(t3, 2)), ')'], ...
       'Location', 'best');
       %'Location', 'best','Interpreter','latex');
title(['Single solution of group GMC at lambda\_ratio = ', num2str(lambda_ratio)])
ax = gca; 
ax.FontSize = 12; 

%%
% figure;
% plot(x1_sg);
% hold on
% plot(x1_gp);
% hold on
% plot(beta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% results for solution path with screening
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the solution path of GMC
% AA2 + FB
t0 = tic;
[xmatrix_sg1, vmatrix_sg1, lambda_seq1] = srls_GMC_path(y, X, 'type', "single",'lambda_min_ratio', 0.01, 'acceleration', "aa2", "screen", false);
t1 = toc(t0);
% AA2 + FB + screning
t0 = tic;
[xmatrix_sg2, vmatrix_sg2, lambda_seq2] = srls_GMC_path(y, X, 'type', "single", 'lambda_min_ratio', 0.01, 'acceleration', "aa2", 'screen', true);
t2 = toc(t0);


norm(xmatrix_sg1-xmatrix_sg2)/norm(xmatrix_sg1)
%% compute the error path
err1 = zeros(1, size(xmatrix_sg1, 1));
err2 = zeros(1, size(xmatrix_sg2, 1));
for i=1: size(xmatrix_sg1, 1)
    err1(i) = norm(xmatrix_sg1(i,:)- beta', 'fro');
    err2(i) = norm(xmatrix_sg2(i,:)- beta', 'fro');
end

%% Plots
figure;
plot(log(lambda_seq1), err1, 'k-', 'LineWidth', 1.5)
hold on; 
plot(log(lambda_seq2), err2, 'b-', 'LineWidth', 1.5)
xlabel('$\log(\lambda)$','Interpreter','latex');
%xlabel('Iteration');
ylabel('Squared error');
legend(['FB + AA2 (t=', num2str(round(t1, 2)), ')'],...
       ['FB + AA2 + screening(t=', num2str(round(t2, 2)), ')'],... 
       ['AA2+FB (t=', num2str(round(t3, 2)), ')'], ...
       'Location', 'best');
       %'Location', 'best','Interpreter','latex');
title(['Error path of GMC'])
ax = gca; 
ax.FontSize = 12; 


%% Compute the solution path of group GMC
% AA2 + FB
t0 = tic;
[xmatrix_gp1, vmatrix_gp1, lambda_seq1] = srls_GMC_path(y, X, 'type', "grouped", 'groups', groups, 'acceleration', "aa2", "screen", false);
t1 = toc(t0);
% AA2 + FB + screning
t0 = tic;
[xmatrix_gp2, vmatrix_gp2, lambda_seq2] = srls_GMC_path(y, X, 'type', "grouped",'groups', groups, 'acceleration', "aa2", 'screen', true);
t2 = toc(t0);


norm(xmatrix_gp1-xmatrix_gp2)/norm(xmatrix_gp1)
%% compute the error path
err1_gp = zeros(1, size(xmatrix_gp1, 1));
err2_gp = zeros(1, size(xmatrix_gp2, 1));
for i=1: size(xmatrix_sg1, 1)
    err1_gp(i) = norm(xmatrix_gp1(i,:)- beta', 'fro');
    err2_gp(i) = norm(xmatrix_gp2(i,:)- beta', 'fro');
end

%% Plots
figure;
plot(log(lambda_seq1), err1_gp, 'k-', 'LineWidth', 1.5)
hold on; 
plot(log(lambda_seq2), err2_gp, 'b-', 'LineWidth', 1.5)
xlabel('$\log(\lambda)$','Interpreter','latex');
%xlabel('Iteration');
ylabel('Squared error');
legend(['FB + AA2 (t=', num2str(round(t1, 2)), ')'],...
       ['FB + AA2 + screening(t=', num2str(round(t2, 2)), ')'],... 
       ['AA2+FB (t=', num2str(round(t3, 2)), ')'], ...
       'Location', 'best');
       %'Location', 'best','Interpreter','latex');
title(['Error path of group GMC'])
ax = gca; 
ax.FontSize = 12; 