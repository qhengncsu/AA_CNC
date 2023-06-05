addpath("Y:\MyDocuments\Xiaoqian\Git\GMC_acceleration_matlab\svt")
%%
%load data_mice
% A = csvread('yeastX.csv',1,1);
% Y = csvread('yeastY.csv',1,1);
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


%% Plot for DY
figure; 
clf;    
plot(log(res_norm_hist_DY), 'b-');
hold on;
plot(log(res_norm_hist_DY_aa), 'r-');
%xlim([-pi pi]);
xlabel('Iteration');
ylabel('Norm of redidual (log scale)');
legend(['Original DY (', num2str(round(t_DY, 2)), ' s)'],...      
       ['AA + DY (',num2str(round(t_DY_aa, 2)), ' s)'], ...
       'FontSize', fsz,...
       'Location', 'best');
       %'Location', 'best','Interpreter','latex');
%title(['DR ', '\lambda = ', num2str(lambda_ratio),'\lambda_{max}'])
% set(gca,'XTick',-3:3); %<- Still need to manually specific tick marks
% set(gca,'YTick',0:10); %<- Still need to manually specific tick marks
ax = gca; 
ax.FontSize = fsz; 

%== save as EPS
%print('portf_DY', '-depsc2', '-r300');

%% Plot for DR
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
%print('portf_DR', '-depsc2', '-r300');



