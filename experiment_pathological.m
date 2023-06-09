params_fixed.splitting = 'FB';
params_fixed.max_iter = 10000;
params_fixed.tol = 1e-5;
params_fixed.early_termination = true;
params_fixed.mem_size = 1;
params_fixed.verbose = true;
params_fixed.printevery = 100;
params_fixed.D = 1e6;
params_fixed.eta = 0;
params_fixed.xi = 0;

[x_opt1, iter1, res_norm_hist1] = fixed_iter(2.1,@forward,@backward,params_fixed,'aa2');

params_fixed.D = 1;
params_fixed.eta = 1e-8;
params_fixed.xi = 1e-14;

[x_opt2, iter2, res_norm_hist2] = fixed_iter(2.1,@forward,@backward,params_fixed,'aa2');
[x_opt3, iter3, res_norm_hist3] = fixed_iter(2.1,@forward,@backward,params_fixed,'original');

plot(1:iter2,res_norm_hist1(1:iter2), 'Color',[0.8,0.58,1], 'LineWidth', 1, 'DisplayName', 'Naive AA')
hold on
plot(1:iter2,res_norm_hist3(1:iter2), 'b-', 'LineWidth', 1, 'DisplayName', 'Original')
hold on
plot(1:iter2,res_norm_hist2(1:iter2), 'r-', 'LineWidth', 1, 'DisplayName', 'Regularized & Safeguarded AA')

l = legend('show','Location','south')
xlabel('iteration');
ylabel('residual norm');
title('Pathological Example')

function x_forward = forward(x)
    x_forward = x - 1/25*grad(x);
end

function x_backward = backward(x)
    x_backward = x;
end

function gradx = grad(x)
    if x<-1
        gradx = x/10-24.9;
    elseif x< 1
        gradx = 25*x;
    else
        gradx = x/10+24.9;
    end
end


