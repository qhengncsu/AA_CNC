rng(2023)

p = 3000;
%Q0 =  gallery('lehmer',p);
Q0 =  gallery('randcorr',p);
%%
gamma = 0.1;
Q = Q0 + gamma* eye(p);

m = randn(p, 1);
r = mean(m);

mu = 0.1;
%%
z0 = ones(3*p, 1);
%%
t0 = tic;
[xhat_aa, res_norm_hist_aa] = cvxmin2(Q, m, r, mu, z0, 'max_iter', 1e4);
t_aa = toc(t0);

%%
t0 = tic;
[xhat, res_norm_hist] = cvxmin2(Q, m, r, mu, z0, 'max_iter', 1e4, 'acceleration', 'original');
t = toc(t0);

%% objective value
Qfun = @(x) Q*x;

z_aa = ( xhat_aa(1:p) +xhat_aa((p+1):(2*p)) + xhat_aa((2*p+1):end))/3;
obj_aa = z_aa'*Qfun(z_aa)/2

z = ( xhat(1:p) +xhat((p+1):(2*p)) + xhat((2*p+1):end))/3;
obj = z'*Qfun(z)/2

err_aa = norm(z_aa-projsplx(z_aa))^2+norm(z_aa-projhalf(z_aa, -m, -r))^2
err = norm(z-projsplx(z))^2+norm(z-projhalf(z, -m, -r))^2

%save portf_DR.mat

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
plot(log(res_norm_hist), 'b-');
hold on;
plot(log(res_norm_hist_aa), 'r-');
%xlim([-pi pi]);
xlabel('Iteration');
ylabel('Norm of redidual (log scale)');
legend(['Original DR (', num2str(round(t, 2)), ' s)'],...      
       ['AA + DR (',num2str(round(t_aa, 2)), ' s)'], ...
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

%% Plot for DY
figure; 
clf;    
plot(log(res_norm_hist_DY), 'b-');
hold on;
plot(log(res_norm_hist_DY_aa), 'r-');
%xlim([-pi pi]);
xlabel('Iteration');
ylabel('Norm of redidual (log scale)');
legend(['Original DY (', num2str(round(t, 2)), ' s)'],...      
       ['AA + DY (',num2str(round(t_aa, 2)), ' s)'], ...
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
