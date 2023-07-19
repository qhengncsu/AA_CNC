function [xhat_matrix, vhat_matrix, intercept, lambda_seq] = srls_GMC_path(y, X, varargin)

% [xhat_matrix, vhat_matrix] = srls_GMC_path(y, X, varargin)
%
% srls_GMC_path: Sparse-Regularized Least Squares with generalized MC (GMC)
% penalty which computes a solution path
%
% Saddle point problem:
%
% argmin_x  argmax_v { F(x,v) =
%  1/2 ||y - A x||^2 + lam ||x||_1 - gamma/2 ||A(x-v)||_2^2 - lam ||v||_1 }
%
% INPUT
%   y 	    response (centered and unit length)
%   X       design matrix (columns centered and unit length)
%
% OUTPUT
%   xhat_matrix, vhat_matrix
% Algorithm: Forward-backward, Theorem 25.8 in Bauschke and Combettes(2011)
% Acceleration: Inertia, Type-II Anderson

params = inputParser;
params.addParameter('type', 'single', @(x) ischar(x)||isstring(x));
params.addParameter('groups', {}, @(x) iscell(x));
params.addParameter('gamma', 0.8, @(x) isnumeric(x));
params.addParameter('splitting', 'FB', @(x) ismember(x,{'DR','FB','FBF'}));
params.addParameter('max_iter', 10000, @(x) isnumeric(x));
params.addParameter('tol_stop', 1e-5, @(x) isnumeric(x));
params.addParameter('lambda_seq', double.empty(0,1), @(x) isvector(x));
params.addParameter('lambda_min_ratio', 0.01, @(x) isnumeric(x));
params.addParameter('screen_off_ratio', 0.05, @(x) isnumeric(x));
params.addParameter('nlambda', 100, @(x) isnumeric(x));
params.addParameter('screen', true, @(x) islogical(x));
params.addParameter('acceleration', 'aa2', @(x) ismember(x,{'original','inertia','aa2'}));
params.addParameter('early_termination', true, @(x) islogical(x));
params.addParameter('mem_size', 5, @(x) isnumeric(x));
params.addParameter('eta', 1e-8, @(x) isnumeric(x));
params.addParameter('printevery', 100, @(x) isnumeric(x));
params.parse(varargin{:});

% soft thresholding
soft = @(x, T) max(1 - T./abs(x), 0) .* x;

type = params.Results.type;
groups = params.Results.groups;
ngroup = length(groups);
gamma = params.Results.gamma;
splitting = params.Results.splitting;
max_iter = params.Results.max_iter;
tol_stop = params.Results.tol_stop;
lambda_min_ratio = params.Results.lambda_min_ratio;
screen_off_ratio = params.Results.screen_off_ratio;
lambda_seq = params.Results.lambda_seq;
nlambda = params.Results.nlambda;
screen = params.Results.screen;
acceleration = params.Results.acceleration;
early_termination = params.Results.early_termination;
mem_size = params.Results.mem_size;
eta = params.Results.eta;
printevery = params.Results.printevery;
params_fixed = struct();
params_fixed.splitting = splitting;
params_fixed.max_iter = max_iter;
params_fixed.tol = tol_stop;
params_fixed.early_termination = early_termination;
params_fixed.mem_size = mem_size;
params_fixed.verbose = false;
params_fixed.eta = eta;
params_fixed.printevery = printevery;
params_fixed.D = 10;
params_fixed.xi = 1e-14;

% data standardization
n = size(X,1);
p = size(X,2);
center = mean(X);
scale = sqrt(sum((X - center).^2)/n);
X = (X - center)./scale;
y = y -mean(y);

Xt = X';
rho = norm(X)^2;
A = @(x) X*x;
AH = @(x) Xt*x;
if strcmp(splitting,'FB')
    mu = 1.99/(rho*(1-2*gamma+2*gamma^2)/(1-gamma));
elseif strcmp(splitting,'FBF')
    gamma_matrix = [1-gamma,gamma;-gamma,gamma];
    mu = 0.99/(rho*norm(gamma_matrix));
end
Xty = Xt*y;
if isempty(lambda_seq)
 
    if strcmp(type,'single')
        lambda_max = max(abs(Xty));
    else
        group_lens = cellfun(@(x) size(x,2), groups);
        Ks = sqrt(group_lens);
        lambda_max = max(group_norm_vec(Xty,groups)./Ks,[],'all');
    end
    lambda_min = lambda_min_ratio*lambda_max;
    inc = -(lambda_max-lambda_min)/(nlambda-1);
    lambda_seq = lambda_max:inc:lambda_min;
  
else
    group_lens = cellfun(@(x) size(x,2), groups);
    Ks = sqrt(group_lens);
    lambda_max = max(lambda_seq);
    lambda_min = min(lambda_seq);

end



% initialization
xhat_matrix = zeros(nlambda,p);
vhat_matrix = zeros(nlambda,p);
intercept = zeros(nlambda, 1);
d_prev = zeros(p,1);
c_prev = -Xty;

for i = 2:nlambda
    lambda = lambda_seq(i);
    % adaptive tolerance?????
    tol_kkt = min(1, lambda/100);  

    xv_current = [xhat_matrix(i-1,:),vhat_matrix(i-1,:)]';
    if screen && (lambda>=screen_off_ratio*lambda_max)
        kkt_fail = true;
        threshold = max(gamma,1-gamma)*(2*lambda - lambda_seq(i-1));
        if strcmp(type,'single')
            keep = ~((abs(d_prev)<threshold) & (abs(c_prev)<threshold));
        else
            keep = ~((group_norm_vec(d_prev,groups)./Ks'<threshold) & (group_norm_vec(c_prev,groups)./Ks'<threshold));
        end
        while kkt_fail
            if strcmp(type,'single')
                a = keep;
            else
                a = false(p,1);
                groups_temp = {};
                l = 0;
                idx = 0;
                for k=1:ngroup
                    if keep(k)
                        a(groups{k}) = true;
                        l = l+1;
                        groups_temp{l} = (idx+1):(idx+length(groups{k}));
                        idx = idx + length(groups{k});
                    end
                end
            end
            p_a = sum(a);
            fprintf('%d features remaining\n', p_a);
            X_a = X(:,a);
            Xt_a = Xt(a,:);
            Xty_a = Xty(a,1);
            A_a = @(x) X_a*x;
            AH_a = @(x) Xt_a*x;
            x_current = xv_current(1:p,1);
            v_current = xv_current((p+1):2*p,1);
            x_current_a = x_current(a,1);
            v_current_a = v_current(a,1);
            xv_current_a = [x_current_a;v_current_a];
            [xv_lambda_a, iter] = fixed_iter(xv_current_a,@forward1,@backward1,params_fixed,acceleration);
            fprintf('lambda = %f solved in %d iterations\n', lambda, iter);
            x_lambda = zeros(p,1);
            v_lambda = zeros(p,1);
            x_lambda(a) = xv_lambda_a(1:p_a);
            v_lambda(a) = xv_lambda_a((p_a+1):2*p_a);
            d_prev = gamma*AH(A(x_lambda-v_lambda));
            c_prev = AH(A(x_lambda)) - Xty - d_prev;
            if strcmp(type,'single')
                kkt_x_case1 = abs(c_prev) > (lambda+tol_kkt);
                kkt_v_case1 = abs(d_prev) > (lambda+tol_kkt);
                kkt_x_case2 = (abs(x_lambda) > tol_kkt) & (abs(c_prev) < lambda-tol_kkt);
                kkt_v_case2 = (abs(v_lambda) > tol_kkt) & (abs(d_prev) < lambda-tol_kkt);
            else
                c_prev_norm = group_norm_vec(c_prev,groups);
                d_prev_norm = group_norm_vec(d_prev,groups);
                x_lambda_norm = group_norm_vec(x_lambda,groups);
                v_lambda_norm = group_norm_vec(v_lambda,groups);
                kkt_x_case1 = c_prev_norm./Ks' > (lambda+tol_kkt);
                kkt_v_case1 = d_prev_norm./Ks' > (lambda+tol_kkt);
                kkt_x_case2 = (x_lambda_norm > tol_kkt) & (c_prev_norm./Ks'<lambda-tol_kkt);
                kkt_v_case2 = (v_lambda_norm > tol_kkt) & (d_prev_norm./Ks'<lambda-tol_kkt);
            end
            kkt_violations = kkt_x_case1 | kkt_v_case1 | kkt_x_case2 | kkt_v_case2;
            kkt_fail = any(kkt_violations);
            total_violations = sum(kkt_violations);
            if kkt_fail
                fprintf('%d optimality conditions violated\n', total_violations);
                keep = keep | kkt_violations;
            end
        end
        
        % unstandardization
        bb = x_lambda;
        xhat_matrix(i,:) = bb./scale';
        intercept(i) = mean(y) - center*bb;
        vhat_matrix(i,:) = v_lambda;
    else
        [xv_lambda, iter] = fixed_iter(xv_current,@forward2,@backward2,params_fixed,acceleration);
        fprintf('lambda = %f solved in %d iterations\n', lambda, iter);
        % unstandardization
        bb = xv_lambda(1:p);
        xhat_matrix(i,:) = bb./scale';
        intercept(i) = mean(y) - center*bb;       
        vhat_matrix(i,:) = xv_lambda((p+1):2*p);
    end
    
 
end

function zxv = forward1(xv)
    x = xv(1:p_a,1);
    v = xv((p_a+1):(2*p_a),1);
    zx = x - mu * ( AH_a(A_a(x + gamma*(v-x))) - Xty_a);
    zv = v - mu * ( gamma * AH_a(A_a(v-x)) );
    zxv = [zx;zv];
end

function xv_next = backward1(zxv)
    zx = zxv(1:p_a,1);
    zv = zxv((p_a+1):(2*p_a),1);
    if strcmp(type,'single')
        x = soft(zx, mu * lambda);
        v = soft(zv, mu * lambda);
    else
        x = soft_group(zx, mu * lambda, groups_temp);
        v = soft_group(zv, mu * lambda, groups_temp);
    end
    xv_next = [x;v];
end

function zxv = forward2(xv)
    x = xv(1:p,1);
    v = xv((p+1):(2*p),1);
    zx = x - mu * ( AH(A(x + gamma*(v-x))) - Xty);
    zv = v - mu * ( gamma * AH(A(v-x)) );
    zxv = [zx;zv];
end

function xv_next = backward2(zxv)
    zx = zxv(1:p,1);
    zv = zxv((p+1):(2*p),1);
    if strcmp(type,'single')
        x = soft(zx, mu * lambda);
        v = soft(zv, mu * lambda);
    else
        x = soft_group(zx, mu * lambda, groups);
        v = soft_group(zv, mu * lambda, groups);
    end
    xv_next = [x;v];
end

function norm_vec = group_norm_vec(x,groups)
    norm_vec = zeros(ngroup,1);
    for j=1:ngroup
        norm_vec(j) = norm(x(groups{j}));
    end
end

function x_prox = soft_group(x,T,groups)
    x_prox = zeros(length(x),1);
    for j=1:length(groups)
        x_prox(groups{j}) = max(1-T*sqrt(length(groups{j}))/norm(x(groups{j})),0)*x(groups{j});
    end
end
end
