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
        while iter<=max_iter
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
        while iter<=max_iter
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
    elseif strcmp(algorithm, 'aa2')
        mem_size = params.mem_size;
        iter = iter + 1;
        mem = x;
        Fmem = F(x);
        while iter<=max_iter
            iter = iter + 1;
            Fmem_mem = Fmem - mem;
            e = ones(size(Fmem_mem,2),1);
            alpha = (Fmem_mem' * Fmem_mem + 1e-10*diag(e)) \ e;
            alp = alpha / sum(alpha);
            x = Fmem * alp;
            Fx = F(x);
            if iter <= mem_size
                mem(:, iter) = x;
                Fmem(:, iter) = Fx;
            else
                mem = [mem(:, 2:end), x];
                Fmem = [Fmem(:, 2:end), Fx];
            end
            res_norm = norm(Fx-x);
            res_norm_hist(iter) = res_norm;
            if mod(iter,100)==0 & verbose
                fprintf('res_norm = %f after %d iterations\n', res_norm, iter);
            end
            if early_termination
                if res_norm/norm(x) < tol
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
    end 
end