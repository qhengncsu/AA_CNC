n = 1e4;
p = 5e3;

%%
rng(2023)
% use sparse matrix
s = 0.1; % sparsity ratio, a alrger sparsity ratio leads to higher speed up ratio
A = randn(n, p);
a = A(:);
a(randsample(n*p, s*n*p)) = 0;
A = reshape(a, [n, p]);

b = randn(n, 1);
mu = 1.99/norm(A'*A);
z0 = lsqr(A, b, 1e-10, 100);
app = 'NNLS';

%%
t0 = tic;
[xhat_aa, res_norm_hist_aa] = cvxmin(A, b, mu, z0, app, 'max_iter', 1e4, 'splitting', 'FB');
t_aa = toc(t0);
%%
t0 = tic;
[xhat, res_norm_hist] = cvxmin(A, b, mu, z0, app, 'max_iter', 1e4, 'acceleration', 'original', 'splitting', 'FB');
t = toc(t0);
%% objective value
obj_aa = norm(A*xhat_aa - b, 'fro')

obj = norm(A*xhat - b, 'fro')

err_aa = sum(min(xhat_aa, 0))
err = sum(min(xhat, 0))

save NNLS_FB.mat

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
plot(log(res_norm_hist), 'k-');
hold on;
plot(log(res_norm_hist_aa), 'r-');
%xlim([-pi pi]);
xlabel('Iteration');
ylabel('Norm of redidual (log scale)');
legend(['Original FB (', num2str(round(t, 2)), ' s)'],...      
       ['AA + FB (',num2str(round(t_aa, 2)), ' s)'], ...
       'FontSize', fsz,...
       'Location', 'best');
       %'Location', 'best','Interpreter','latex');
%title(['DR ', '\lambda = ', num2str(lambda_ratio),'\lambda_{max}'])
% set(gca,'XTick',-3:3); %<- Still need to manually specific tick marks
% set(gca,'YTick',0:10); %<- Still need to manually specific tick marks
ax = gca; 
ax.FontSize = fsz; 

%== save as EPS
print('NNLS_FB', '-depsc2', '-r300');
