function [x, iter, res_norm_hist] = fixed_iter(x0, F, params, algorithm)
    max_iter = params.max_iter;
    verbose = params.verbose;
    tol = params.tol;
    cg_check = 0;
    res_norm_hist = zeros(max_iter,1);
    early_termination = params.early_termination;
    x = x0;
    iter = 0;
    if strcmp(algorithm, 'original')
        while iter<=max_iter-1
            iter = iter + 1;
            Fx = F(x);
            res = Fx - x;
            res_norm = norm(res);
            res_norm_hist(iter) = res_norm;
            if mod(iter,100)==0 & verbose
                fprintf('res_norm = %f after %d iterations\n', res_norm, iter);
            end
            x = Fx;
            if early_termination
                if res_norm/norm(x) < tol
                    cg_check = cg_check + 1;
                else
                    cg_check = 0;
                end
                if cg_check >= 10
                    break
                end
            end
        end     
    elseif strcmp(algorithm, 'nesterov')
        k = 1;
        x_prev = x;
        x_prev2 = x;
        res_norm_prev = Inf;
        while iter<=max_iter-1
            iter = iter + 1;
            x_tilde = x_prev + (k-1)/(k+2)*(x_prev - x_prev2);
            Fx = F(x_tilde);
            res_norm = norm(Fx-x_tilde);
            if mod(iter,100)==0 & verbose
                fprintf('res_norm = %f after %d iterations\n', res_norm, iter);
            end
            x_prev2 = x_prev;
            x_prev = Fx;
            k = k + 1;
            if res_norm > res_norm_prev
                k = 1;
                x_prev = x_prev2;
                res_norm_hist(iter) = res_norm_prev;
            else
                res_norm_prev = res_norm;
                res_norm_hist(iter) = res_norm;
            end
            if early_termination
                if res_norm/norm(x_tilde) < tol
                    cg_check = cg_check + 1;
                else
                    cg_check = 0;
                end
                if cg_check >= 10
                    x = Fx;
                    break
                end
            end
        end
    elseif strcmp(algorithm, 'inertia')
        x_prev = x;
        x_prev2 = x;
        while iter<=max_iter-1
            iter = iter + 1;
            x_tilde = x_prev + 0.333*(x_prev - x_prev2);
            Fx = F(x_tilde);
            res_norm = norm(Fx-x_tilde);
            if mod(iter,100)==0 & verbose
                fprintf('res_norm = %f after %d iterations\n', res_norm, iter);
            end
            x_prev2 = x_prev;
            x_prev = Fx;
            res_norm_hist(iter) = res_norm;
            if early_termination
                if res_norm/norm(x_tilde) < tol
                    cg_check = cg_check + 1;
                else
                    cg_check = 0;
                end
                if cg_check >= 10
                    x = Fx;
                    break
                end
            end
        end
    elseif strcmp(algorithm, 'aa2')
        mem_size = params.mem_size;
        eta = params.eta;
        xk_1 = x0;
        xk = F(xk_1);
        Xk = xk;
        gk_1 = xk_1 - xk;
        while iter<=max_iter-1
            iter = iter + 1;
            xkp1_FB = F(xk);
            gk = xk - xkp1_FB;
            %populate the history before starting anderson iterations
            if iter <= mem_size    
                Xk(:, iter+1) = xkp1_FB;
                Yk(:, iter) = gk - gk_1;
                Sk(:, iter) = xk - xk_1;
                xk_1 = xk;
                xk = xkp1_FB;
                gk_1 = gk;
            else
                Yk = [Yk(:, 2:end), gk - gk_1];
                Sk = [Sk(:, 2:end), xk - xk_1];
                Xk = [Xk(:, 2:end), xkp1_FB];
                e = ones(size(Yk,2),1);
                Yk_norm2 = norm(Yk,"fro")^2;
                Sk_norm2 = norm(Sk,"fro")^2;
                gamma_k = (Yk' * Yk + (eta*(Sk_norm2+Yk_norm2)+1e-14)*diag(e))\(Yk' * gk);
                alpha_k = [gamma_k(1);gamma_k(2:end)-gamma_k(1:(end-1));1-gamma_k(end)];
                xkp1_aa = Xk*alpha_k;
                xk = xkp1_aa;
            end
            xk_1 = xk;
            gk_1 = gk;
            res_norm = norm(gk_1);
            res_norm_hist(iter) = res_norm;
            if mod(iter,100)==0 & verbose
                fprintf('res_norm = %f after %d iterations\n', res_norm, iter);
            end
            if early_termination
                if res_norm/norm(xk_1) < tol
                    cg_check = cg_check + 1;
                else
                    cg_check = 0;
                end
                if cg_check >= 10
                    x = xk;
                    break
                end
            end
        end
        x = xk;
    end 
end