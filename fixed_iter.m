function [z, iter, res_norm_hist] = fixed_iter(z0, forward, backward, params, acceleration)
    max_iter = params.max_iter;
    verbose = params.verbose;
    splitting = params.splitting;
    tol = params.tol;
    cg_check = 0;
    res_norm_hist = zeros(max_iter,1);
    early_termination = params.early_termination;
    z_length = length(z0);
    printevery = params.printevery;
    if strcmp(splitting,'DY')
        projection = params.projection;
    elseif strcmp(splitting,'DK')
        u = -z0;
    end
    z = z0;
    iter = 0;
    if strcmp(acceleration, 'original')
        while iter<=max_iter-1
            iter = iter + 1;
            if strcmp(splitting, 'FB')
                Fz = backward(forward(z));
                res = z - Fz;
            elseif strcmp(splitting, 'FBF')
                z_f = forward(z);
                z_fb = backward(z_f);
                z_fbf = forward(z_fb);
                res = z_f - z_fbf;
                Fz = z - res;
            elseif strcmp(splitting,'DY')
                z_Q = backward(z);
                z_R = projection(2.*z_Q - z- forward(z_Q));
                Fz = z-z_Q+z_R;
                res = z - Fz;
            elseif strcmp(splitting,'DK')
                z_Q = z + backward(-z);
                Fz = z_Q - forward(z_Q+u);
                res = z - Fz;
            end
            res_norm = norm(res);
            res_norm_hist(iter) = res_norm;
            if mod(iter,printevery)==0 && verbose
                fprintf('res_norm = %f after %d iterations\n', res_norm, iter);
            end
            z = Fz;
            if early_termination
                if res_norm/(norm(z)+1) < tol
                    cg_check = cg_check + 1;
                else
                    cg_check = 0;
                end
                if cg_check >= 10
                    break
                end
            end
        end     
    % elseif strcmp(acceleration, 'nesterov')
    %     k = 1;
    %     z_prev = z;
    %     z_prev2 = z;
    %     res_norm_prev = Inf;
    %     while iter<=max_iter-1
    %         iter = iter + 1;
    %         z_tilde = z_prev + (k-1)/(k+2)*(z_prev - z_prev2);
    %         Fz = backward(forward(z_tilde));
    %         res_norm = norm(Fz-z_tilde);
    %         if mod(iter,100)==0 && verbose
    %             fprintf('res_norm = %f after %d iterations\n', res_norm, iter);
    %         end
    %         z_prev2 = z_prev;
    %         z_prev = Fz;
    %         k = k + 1;
    %         if res_norm > res_norm_prev
    %             k = 1;
    %             z_prev = z_prev2;
    %             res_norm_hist(iter) = res_norm_prev;
    %         else
    %             res_norm_prev = res_norm;
    %             res_norm_hist(iter) = res_norm;
    %         end
    %         if early_termination
    %             if res_norm/norm(z_tilde) < tol
    %                 cg_check = cg_check + 1;
    %             else
    %                 cg_check = 0;
    %             end
    %             if cg_check >= 5
    %                 z = Fz;
    %                 break
    %             end
    %         end
    %     end
    elseif strcmp(acceleration, 'inertia')
        z_prev = z;
        z_prev2 = z;
        while iter<=max_iter-1
            iter = iter + 1;
            z_tilde = z_prev + 0.333*(z_prev - z_prev2);
            if strcmp(splitting, 'FB')
                Fz = backward(forward(z_tilde));
            else
                error("Not Implemented.");
            end
            res_norm = norm(Fz-z_tilde);
            if mod(iter,printevery)==0 && verbose
                fprintf('res_norm = %f after %d iterations\n', res_norm, iter);
            end
            z_prev2 = z_prev;
            z_prev = Fz;
            res_norm_hist(iter) = res_norm;
            if early_termination
                if res_norm/norm(z_tilde) < tol
                    cg_check = cg_check + 1;
                else
                    cg_check = 0;
                end
                if cg_check >= 10
                    z = Fz;
                    break
                end
            end
        end
    elseif strcmp(acceleration, 'aa2')
        mem_size = params.mem_size;
        eta = params.eta;
        zk_1 = z0;
        if strcmp(splitting, 'FB')
            zk = backward(forward(zk_1));
        elseif strcmp(splitting, 'FBF')
            z_f = forward(zk_1);
            z_fb = backward(z_f);
            z_fbf = forward(z_fb);
            zk = zk_1 - z_f + z_fbf; 
        elseif strcmp(splitting, 'DY')
            z_Q = backward(zk_1);
            z_R = projection(2.*z_Q - zk_1- forward(z_Q));
            zk = zk_1-z_Q+z_R;
        elseif strcmp(splitting, 'DK')
            z_Q = zk_1 + backward(-zk_1);
            zk = z_Q - forward(z_Q+u);
        end
        Zk(:,1) = zk;
        gk_1 = zk_1 - zk;
        res_norm = norm(gk_1);
        res_norm_hist(1) = res_norm;
        g0 = gk_1;
        norm_g0 = norm(g0);
        D = 1e8;
        epsilon = 1e-8;
        i = 0;
        while iter<=max_iter-2
            iter = iter + 1;
            if strcmp(splitting, 'FB')
                zkp1_OS = backward(forward(zk));
                gk = zk - zkp1_OS;
            elseif strcmp(splitting, 'FBF')
                z_f = forward(zk);
                z_fb = backward(z_f);
                z_fbf = forward(z_fb);
                gk_fb = zk - z_fb;
                gk = z_f - z_fbf;
                zkp1_OS = zk - gk;
            elseif strcmp(splitting,'DY')
                z_Q = backward(zk);
                z_R = projection(2.*z_Q - zk- forward(z_Q));
                zkp1_OS = zk-z_Q+z_R;
                gk = zk-zkp1_OS;
            elseif strcmp(splitting,'DK')
                z_Q = zk + backward(-zk);
                zkp1_OS = z_Q - forward(z_Q+u);
                gk = zk-zkp1_OS;
            end
            if iter <= mem_size    
                Zk(:, iter+1) = zkp1_OS;
                Yk(:, iter) = gk - gk_1;
                Sk_norms(iter) = norm(zk - zk_1);
            else
                Yk = [Yk(:, 2:end), gk - gk_1];
                Sk_norms = [Sk_norms(:, 2:end), norm(zk - zk_1)];
                Zk = [Zk(:, 2:end), zkp1_OS];
            end
            Yk_norm2 = norm(Yk,"fro")^2;
            Sk_norm2 = sum(Sk_norms,'all');
            left = Yk' * Yk + (eta*(Sk_norm2+Yk_norm2)+1e-14)*eye(min(mem_size,iter));
            gamma_k = left\(Yk' * gk);
            alpha_k = [gamma_k(1);gamma_k(2:end)-gamma_k(1:(end-1));1-gamma_k(end)];
            zkp1_aa = Zk*alpha_k;
            if ismember(splitting,{'FB','DR','DY','DK'}) && norm(gk)<=D*norm_g0*(i+1)^(-1-epsilon)
                zk = zkp1_aa;
                i = i+1;
            elseif strcmp(splitting, 'FBF') && norm(gk_fb)<=0.5*D*norm_g0*(i+1)^(-1-epsilon)
                zk = zkp1_aa;
                i = i+1;
            else
                zk = zkp1_OS;
            end
            zk_1 = zk;
            gk_1 = gk;
            res_norm = norm(gk_1);
            res_norm_hist(iter+1) = res_norm;
            if mod(iter,printevery)==0 && verbose
                fprintf('res_norm = %f after %d iterations\n', res_norm, iter);
            end
            if early_termination
                if res_norm/(norm(zk_1)+1) < tol
                    cg_check = cg_check + 1;
                else
                    cg_check = 0;
                end
                if cg_check >= 10
                    break
                end
            end
        end
        z = zk;
    end
    if strcmp(splitting, 'DY')
        z = backward(z);
    elseif strcmp(splitting, 'FBF')
        z = backward(forward(z));
    elseif strcmp(splitting, 'DK')
        z = backward(-z);
    end
end
