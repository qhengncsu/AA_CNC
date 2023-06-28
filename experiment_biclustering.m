U = csvread("data/lung500.csv",1,1);
times_original = zeros(5,1);
times_aa2 = zeros(5,1);
lambdas = 4:0.5:6;
for i=1:5
    lambda = 10^(lambdas(i));
    tic
        if i==2
            [xhat1, res_norm_hist1] = dykstra_splitting(U, 'biclustering', lambda, 'acceleration', 'original','printevery',1);
        else 
            [xhat1, ~] = dykstra_splitting(U, 'biclustering', lambda, 'acceleration', 'original','printevery',1);
        end   
    times_original(i) = toc;
    tic
        if i==2
            [xhat2, res_norm_hist2] = dykstra_splitting(U, 'biclustering', lambda, 'acceleration', 'aa2','printevery',1);
        else 
            [xhat2, ~] = dykstra_splitting(U, 'biclustering', lambda, 'acceleration', 'aa2','printevery',1);
        end 
    times_aa2(i) = toc;
end
save("results/biclustering.mat","times_original","times_aa2","res_norm_hist1","res_norm_hist2")
