%% add path to the code 
%addpath("Y:\MyDocuments\Xiaoqian\GMC-computation\code")
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

groups = cell(gp_num,1);
for i=1:gp_num
    groups{i} = ((i-1)*gp_size+1):(i*gp_size);
end

% on a single lambda value, report the iterations
lambda_ratio = 0.1;
gamma = 0.8;
%% algorithm comparison on group GMC: FB v.s. FB+AA

% vanilla
t0 = tic;
[x_fb, v_fb, res_norm_fb] = srls_GMC_acc(y, X, lambda_ratio, 'type', 'group','groups',groups,...
                  'acceleration', 'original',"gamma",gamma, 'splitting', 'FB');
t1 = toc(t0);


% AA
t0 = tic;
[x_fbaa, v_fbaa, res_norm_fbaa] = srls_GMC_acc(y, X,  lambda_ratio, 'type','group','groups',groups,...
                  'acceleration', "aa2","gamma",gamma, 'splitting', 'FB');
t2 = toc(t0);
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
plot(log(res_norm_fb), 'k-');
hold on;
plot(log(res_norm_fbaa), 'r-');
%xlim([-pi pi]);
xlabel('Iteration');
ylabel('Norm of redidual (log scale)');
legend(['Original FB (', num2str(round(t1, 2)), ' s)'],...      
       ['AA + FB (',num2str(round(t2, 2)), ' s)'], ...
       'FontSize', fsz,...
       'Location', 'best');
       %'Location', 'best','Interpreter','latex');
title(['Group GMC, ', '\lambda = ', num2str(lambda_ratio),'\lambda_{max}'])
% set(gca,'XTick',-3:3); %<- Still need to manually specific tick marks
% set(gca,'YTick',0:10); %<- Still need to manually specific tick marks
ax = gca; 
ax.FontSize = fsz; 

%== save as EPS
print('cvg_group_FB', '-depsc2', '-r300');

%%
%% algorithm comparison on group GMC: FBF v.s. FBF+AA
% vanilla
t0 = tic;
[x_fbf, v_fbf, res_norm_fbf] = srls_GMC_acc(y, X, lambda_ratio, 'type', 'group','groups',groups,...
                  'acceleration', 'original',"gamma",gamma, 'splitting', 'FBF');
t3 = toc(t0);

% AA
t0 = tic;
[x_fbfaa, v_fbfaa, res_norm_fbfaa] = srls_GMC_acc(y, X,  lambda_ratio, 'type','group','groups',groups,...
                  'acceleration', "aa2","gamma",gamma, 'splitting', 'FBF');
t4 = toc(t0);
%% Plot for FBF
figure; 
clf;    
plot(log(res_norm_fbf), 'k-');
hold on;
plot(log(res_norm_fbfaa), 'r-');
%xlim([-pi pi]);
xlabel('Iteration');
ylabel('Norm of redidual (log scale)');
legend(['Original FBF (', num2str(round(t3, 2)), ' s)'],...      
       ['AA + FBF (',num2str(round(t4, 2)), ' s)'], ...
       'FontSize', fsz,...
       'Location', 'best');
       %'Location', 'best','Interpreter','latex');
title(['Group GMC, ', '\lambda = ', num2str(lambda_ratio),'\lambda_{max}'])
% set(gca,'XTick',-3:3); %<- Still need to manually specific tick marks
% set(gca,'YTick',0:10); %<- Still need to manually specific tick marks
ax = gca; 
ax.FontSize = fsz; 

%== save as EPS
print('cvg_group_FBF', '-depsc2', '-r300');



%%
%% algorithm comparison on GMC: DR v.s. DR+AA
% vanilla
t0 = tic;
[x_dr, v_dr, res_norm_dr] = srls_GMC_acc(y, X, lambda_ratio, 'type', 'group','groups',groups,...
                  'acceleration', 'original',"gamma",gamma, 'splitting', 'DR');
t5 = toc(t0);

% AA
t0 = tic;
[x_draa, v_draa, res_norm_draa] = srls_GMC_acc(y, X, lambda_ratio, 'type','group', 'groups', groups,...
                  'acceleration', "aa2","gamma",gamma, 'splitting', 'DR');
t6 = toc(t0);
%% Plot for DR
figure; 
clf;    
plot(log(res_norm_dr), 'k-');
hold on;
plot(log(res_norm_draa), 'r-');
%xlim([-pi pi]);
xlabel('Iteration');
ylabel('Norm of redidual (log scale)');
legend(['Original DR (', num2str(round(t5, 2)), ' s)'],...      
       ['AA + DR (',num2str(round(t6, 2)), ' s)'], ...
       'FontSize', fsz,...
       'Location', 'best');
       %'Location', 'best','Interpreter','latex');
title(['Group GMC, ', '\lambda = ', num2str(lambda_ratio),'\lambda_{max}'])
% set(gca,'XTick',-3:3); %<- Still need to manually specific tick marks
% set(gca,'YTick',0:10); %<- Still need to manually specific tick marks
ax = gca; 
ax.FontSize = fsz; 

%== save as EPS
print('cvg_group_DR', '-depsc2', '-r300');

%%

save cvg_group.mat