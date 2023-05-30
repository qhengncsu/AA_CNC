%% add path to the code 
addpath("Y:\MyDocuments\Xiaoqian\GMC-computation\code")
%% data generation
n = 1000;
p = 10000;
SNR = 1;
rng(2023)
rho = 0;
Sigma = zeros(p, p);
for k=1:p
    for l=1:p
        Sigma(k, l) = rho^(abs(k-l));
    end
end
X = mvnrnd(zeros(p,1),Sigma, n);

beta = [ones(p/100,1);  zeros(p*99/100,1)];
y = X*beta + sqrt(beta'*Sigma*beta/SNR)*randn(n,1);


%% algorithm comparison: vanilla, Nesterov, AA on GMC
% on a single lambda value, report the iterations
lambda_ratio = 0.1;
gamma = 0.8;
% vanilla
t0 = tic;
[x1_sg, v1_sg, res_norm_sg1] = srls_GMC_acc(y, X, lambda_ratio, 'type', 'single',...
                  'acceleration', 'original',"gamma",gamma, 'splitting', 'FB');
t1 = toc(t0);

% % Nesterov
% t0 = tic;
% [x2_sg, v2_sg, res_norm_sg2] = srls_GMC_acc(y, X, lambda_ratio, 'type', "single", 'acceleration', "inertia","gamma",gamma);
% t2 = toc(t0);

% AA
t0 = tic;
[x3_sg, v3_sg, res_norm_sg3] = srls_GMC_acc(y, X,  lambda_ratio, 'type','single',...
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
plot(log(res_norm_sg1), 'k-');
hold on;
plot(log(res_norm_sg3), 'r-');
%xlim([-pi pi]);
xlabel('Iteration');
ylabel('Norm of redidual (log scale)');
legend(['Original FB (', num2str(round(t1, 2)), ' s)'],...      
       ['AA + FB (',num2str(round(t2, 2)), ' s)'], ...
       'FontSize', fsz,...
       'Location', 'best');
       %'Location', 'best','Interpreter','latex');
title(['FB ', '\lambda = ', num2str(lambda_ratio),'\lambda_{max}'])
% set(gca,'XTick',-3:3); %<- Still need to manually specific tick marks
% set(gca,'YTick',0:10); %<- Still need to manually specific tick marks
ax = gca; 
ax.FontSize = fsz; 

%== save as EPS
print('cvg_single_FB', '-depsc2', '-r300');

%%
%% algorithm comparison: vanilla, Nesterov, AA on GMC
% vanilla
t0 = tic;
[x1_sg, v1_sg, res_norm_sg1] = srls_GMC_acc(y, X, lambda_ratio, 'type', 'single',...
                  'acceleration', 'original',"gamma",gamma, 'splitting', 'FBF');
t3 = toc(t0);

% AA
t0 = tic;
[x3_sg, v3_sg, res_norm_sg3] = srls_GMC_acc(y, X,  lambda_ratio, 'type','single',...
                  'acceleration', "aa2","gamma",gamma, 'splitting', 'FBF');
t4 = toc(t0);
%% Plot for FBF
figure; 
clf;    
plot(log(res_norm_sg1), 'k-');
hold on;
plot(log(res_norm_sg3), 'r-');
%xlim([-pi pi]);
xlabel('Iteration');
ylabel('Norm of redidual (log scale)');
legend(['Original FBF (', num2str(round(t3, 2)), ' s)'],...      
       ['AA + FBF (',num2str(round(t4, 2)), ' s)'], ...
       'FontSize', fsz,...
       'Location', 'best');
       %'Location', 'best','Interpreter','latex');
title(['FBF, ', '\lambda = ', num2str(lambda_ratio),'\lambda_{max}'])
% set(gca,'XTick',-3:3); %<- Still need to manually specific tick marks
% set(gca,'YTick',0:10); %<- Still need to manually specific tick marks
ax = gca; 
ax.FontSize = fsz; 

%== save as EPS
print('cvg_single_FBF', '-depsc2', '-r300');