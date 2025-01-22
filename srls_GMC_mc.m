function [xhat, vhat, res_norm_hist] = srls_GMC_mc(y, lambda, varargin)
% [xhat, vhat, res_norm_hist] = srls_GMC_matrix(y, app, lambda, varargin)
% INPUT
%   y 	    observed matrix
%   lambda  penalty parameter
% OUTPUT
%   xhat, vhat, res_norm_hist
params = inputParser;
%params.addParameter('H', fspecial('average',3), @(x) isnumeric(x));
params.addParameter('gamma', 0.8, @(x) isnumeric(x));
params.addParameter('splitting', 'FBF', @(x) ismember(x,{'FB','FBF'}));
params.addParameter('max_iter', 10000, @(x) isnumeric(x));
params.addParameter('tol_stop', 1e-5, @(x) isnumeric(x));
params.addParameter('early_termination', true, @(x) islogical(x));
params.addParameter('acceleration', 'aa2', @(x) ismember(x,{'original','aa2'}));
params.addParameter('mask', ones(size(y)), @(x) isnumeric(x));
params.addParameter('mem_size', 10, @(x) isnumeric(x));
params.addParameter('eta', 1e-2, @(x) isnumeric(x));
params.addParameter('printevery', 100, @(x) isnumeric(x));
%params.addParameter('lower', 0, @(x) isnumeric(x));
%params.addParameter('upper', 1, @(x) isnumeric(x));
params.addParameter('xv0', [full(y(:));0*full(y(:))], @(x) isnumeric(x));
params.parse(varargin{:});
%H = params.Results.H;
gamma = params.Results.gamma;
splitting = params.Results.splitting;
max_iter = params.Results.max_iter;
tol_stop = params.Results.tol_stop;
early_termination = params.Results.early_termination;
acceleration = params.Results.acceleration;
mask = params.Results.mask;
mem_size = params.Results.mem_size;
eta = params.Results.eta;
printevery = params.Results.printevery;
%lower = params.Results.lower;
%upper = params.Results.upper;
xv0 = params.Results.xv0;
params_fixed = struct();
params_fixed.splitting = splitting;
params_fixed.max_iter = max_iter;
params_fixed.tol = tol_stop;
params_fixed.early_termination = early_termination;
params_fixed.mem_size = mem_size;
params_fixed.verbose = true;
params_fixed.eta = eta;
params_fixed.printevery = printevery;
params_fixed.D = 10;
params_fixed.xi = 1e-14;
n1 = size(y,1);
n2 = size(y,2);
rho = 1;
n_total = n1*n2;

%if strcmp(app,'deblurring')
    %A = @(x) imfilter(x,H,'circular');
    %AH = @(x) imfilter(x,H,'circular');
    %AHy = AH(y);
    %AHA = @(x) AH(A(x));
%elseif strcmp(app,'matrix completion')
AHA = @(x) mask.*x;
AHy = AHA(y);
    %projection = @(x) [min(max(x(1:n_total),lower),upper);x((n_total+1):end)];
    %params_fixed.projection = projection;
%end

if strcmp(splitting,'FBF') 
    gamma_matrix = [1-gamma,gamma;-gamma,gamma];
    mu = 0.99/(rho*norm(gamma_matrix))
else
    mu = 1.99*min(1,(1-gamma)/gamma)/rho;
end
[xv_lambda, iter, res_norm_hist] = fixed_iter(xv0,@forward,@backward,params_fixed,acceleration);
xhat = reshape(xv_lambda(1:n_total),[n1,n2]);
vhat = reshape(xv_lambda((n_total+1):(2*n_total)), [n1,n2]);
fprintf('lambda = %f solved in %d iterations\n', lambda, iter);

function zxv = forward(xv)
    x = reshape(xv(1:n_total),[n1,n2]);
    v = reshape(xv((n_total+1):(2*n_total)), [n1,n2]);
    zx = x - mu * (AHA(x+gamma*(v-x))-AHy);
    zv = v - mu * (gamma*AHA(v-x));
    zxv = [zx(:);zv(:)];
end

function xv_next = backward(zxv)
    % if strcmp(app,'matrix completion')
        % if gamma>0
    zx = reshape(zxv(1:n_total),[n1,n2]);
    zv = reshape(zxv((n_total+1):(2*n_total)), [n1,n2]);
    x = svt(zx, mu * lambda);
	v = svt(zv, mu * lambda);
    xv_next = [x(:);v(:)];
%         else
%             zx = reshape(zxv,[n1,n2]);
%             x = svt(zx, mu * lambda);
%             xv_next = x(:);
%         end
%     elseif strcmp(app,'deblurring')
%         if gamma>0
%             zx = reshape(zxv(1:n_total),[n1,n2]);
%             zv = reshape(zxv((n_total+1):(2*n_total)), [n1,n2]);
%             x = chambolle_prox_TV_stop(zx,'lambda',mu * lambda,'maxiter',50);
%         	v = chambolle_prox_TV_stop(zv,'lambda',mu * lambda,'maxiter',50);
%             xv_next = [x(:);v(:)];
%         else
%             zx = reshape(zxv,[n1,n2]);
%             x = chambolle_prox_TV_stop(zx,'lambda',mu * lambda,'maxiter',50);
%             xv_next = x(:);
%         end
%     end
% end
end
end