%% add path to the code 
addpath('/rsrch6/home/bcb/xliu31/GMC/code')
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

%% Plot settings
width = 4;     % Width in inches 
height = 3;    % Height in inches 
alw = 2;    % AxesLineWidth 
fsz = 11;      % Fontsize 
lw = 2;      % LineWidth 
msz = 2;       % MarkerSize 

% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]); 

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);

%% Plot for FB
figure; 
clf;    
plot(log(res_norm_hist_DR), 'b-');
hold on;
plot(log(res_norm_hist_DR_aa), 'r-');
%xlim([-pi pi]);
xlabel('Iteration');
ylabel('Norm of redidual (log scale)');
legend(['Original DR (', num2str(round(t_DR, 2)), ' s)'],...      
       ['AA + DR (',num2str(round(t_DR_aa, 2)), ' s)'], ...
       'FontSize', fsz,...
       'Location', 'best');
       %'Location', 'best','Interpreter','latex');
%title(['DR ', '\lambda = ', num2str(lambda_ratio),'\lambda_{max}'])
% set(gca,'XTick',-3:3); %<- Still need to manually specific tick marks
% set(gca,'YTick',0:10); %<- Still need to manually specific tick marks
ax = gca; 
ax.FontSize = fsz; 

%== save as EPS
print('NNLS_DR_n3e4_s5', '-depsc2', '-r300');

%% Plot for FB
figure; 
clf;    
plot(log(res_norm_hist_FB), 'b-');
hold on;
plot(log(res_norm_hist_FB_aa), 'r-');
%xlim([-pi pi]);
xlabel('Iteration');
ylabel('Norm of redidual (log scale)');
legend(['Original FB (', num2str(round(t_FB, 2)), ' s)'],...      
       ['AA + FB (',num2str(round(t_FB_aa, 2)), ' s)'], ...
       'FontSize', fsz,...
       'Location', 'best');
       %'Location', 'best','Interpreter','latex');
%title(['DR ', '\lambda = ', num2str(lambda_ratio),'\lambda_{max}'])
% set(gca,'XTick',-3:3); %<- Still need to manually specific tick marks
% set(gca,'YTick',0:10); %<- Still need to manually specific tick marks
ax = gca; 
ax.FontSize = fsz; 

%== save as EPS
print('NNLS_FB_n3e4_s5', '-depsc2', '-r300');


%%
 save NNLS_results_n3e4_s5.mat
