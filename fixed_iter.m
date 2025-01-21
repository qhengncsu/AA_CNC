function [zstar, iter, res_norm_hist] = fixed_iter(z0, forward, backward, params, acceleration)
    max_iter = params.max_iter;
    verbose = params.verbose;
    splitting = params.splitting;
    tol = params.tol;
    res_norm_hist = zeros(max_iter,1);
    early_termination = params.early_termination;
    printevery = params.printevery;
    if strcmp(splitting,'DY')
        projection = params.projection;
    end
    z = z0;
    iter = 0;
    if strcmp(acceleration, 'original')
        while iter<=max_iter-1
            iter = iter + 1;
            if strcmp(splitting, 'FB')
                Fz = backward(forward(z));
            elseif strcmp(splitting, 'FBF')
                z_f = forward(z);
                z_fb = backward(z_f);
                z_fbf = forward(z_fb);
                Fz = z - z_f + z_fbf;
            elseif strcmp(splitting,'DY')
                z_R = projection(z);
                z_Q = backward(z_R - z + forward(z_R));
                Fz = z-z_R+z_Q;
            elseif strcmp(splitting,'DR')
                z_P = forward(z);
                z_Q = backward(2*z_P-z);
                Fz = z - z_P + z_Q;
            end
            res = z - Fz;
            res_norm = norm(res);
            res_norm_hist(iter) = res_norm;
            if mod(iter,printevery)==0 && verbose
                fprintf('res_norm = %f after %d iterations\n', res_norm, iter);
            end
            z = Fz;
            if early_termination
                if res_norm/(norm(z)+1) < tol
                    break
                end
            end
        end
        zstar = z;
    elseif strcmp(acceleration, 'aa2')
        mem_size = params.mem_size;
        eta = params.eta;
        D = params.D;
        xi = params.xi;
        zk_1 = z0;
        if strcmp(splitting, 'FB')
            zk = backward(forward(zk_1));
        elseif strcmp(splitting, 'FBF')
            z_f = forward(zk_1);
            z_fb = backward(z_f);
            z_fbf = forward(z_fb);
            zk = zk_1 - z_f + z_fbf; 
        elseif strcmp(splitting, 'DY')
            z_R = projection(zk_1);
            z_Q = backward(z_R - zk_1 + forward(z_R));
            zk = zk_1-z_R+z_Q;
        elseif strcmp(splitting,'DR')
            z_P = forward(zk_1);
            z_Q = backward(2*z_P-zk_1);
            zk = zk_1 - z_P + z_Q;
        end
        Zk(:,1) = zk;
        gk_1 = zk_1 - zk;
        res_norm = norm(gk_1);
        res_norm_hist(1) = res_norm;
        g0 = gk_1;
        norm_g0 = norm(g0);
        epsilon = 1e-6;
        i = 0;
        total_safeguards = 0;
        while iter<=max_iter-1
            iter = iter + 1;
            if strcmp(splitting, 'FB')
                zkp1_OS = backward(forward(zk));
            elseif strcmp(splitting, 'FBF')
                z_f = forward(zk);
                z_fb = backward(z_f);
                z_fbf = forward(z_fb);
                gk_fb = zk - z_fb;
                zkp1_OS = zk - z_f + z_fbf;
            elseif strcmp(splitting,'DY')
                z_R = projection(zk);
                z_Q = backward(z_R - zk + forward(z_R));  
                zkp1_OS = zk-z_R+z_Q;
            elseif strcmp(splitting,'DR')
                z_P = forward(zk);
                z_Q = backward(2*z_P-zk);
                zkp1_OS = zk - z_P + z_Q;
            end
            gk = zk - zkp1_OS;
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
            left = Yk' * Yk + (eta*(Sk_norm2+Yk_norm2)+xi)*eye(min(mem_size,iter));
            gamma_k = left\(Yk' * gk);
            alpha_k = [gamma_k(1);gamma_k(2:end)-gamma_k(1:(end-1));1-gamma_k(end)];
            zkp1_aa = Zk*alpha_k;
            if ismember(splitting,{'FB','DR','DY'}) && norm(gk)<=D*norm_g0*(i+1)^(-1-epsilon)
                zk = zkp1_aa;
                i = i+1;
            elseif strcmp(splitting, 'FBF') && norm(gk_fb)<=0.5*D*norm_g0*(i+1)^(-1-epsilon)
                zk = zkp1_aa;
                i = i+1;
            else
                zk = zkp1_OS;
                total_safeguards = total_safeguards + 1;
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
                    break
                end
            end
        end
        zstar = zk;
        fprintf('Safeguard invoked %d times!\n', total_safeguards);
    end
    if strcmp(splitting, 'DY')
        zstar = projection(zstar);
    elseif strcmp(splitting, 'DR')
        zstar = forward(zstar);
    end
end
