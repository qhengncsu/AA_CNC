U = csvread("lung500.csv",1,1);
times_original = zeros(5,1);
times_aa2 = zeros(5,1);
lambdas = 4:0.5:6;
for i=1:5
    lambda = 10^(lambdas(i));
    tic
        if i==3
            [xhat1, res_norm_hist1] = dykstra_splitting(U, 'biclustering', lambda, 'acceleration', 'original','printevery',1);
        else 
            [xhat1, ~] = dykstra_splitting(U, 'biclustering', lambda, 'acceleration', 'original','printevery',1);
        end   
    times_original(i) = toc;
    tic
        if i==3
            [xhat2, res_norm_hist2] = dykstra_splitting(U, 'biclustering', lambda, 'acceleration', 'aa2','printevery',1);
        else 
            [xhat2, ~] = dykstra_splitting(U, 'biclustering', lambda, 'acceleration', 'aa2','printevery',1);
        end 
    times_aa2(i) = toc;
end
csvwrite('times_biclustering_original.csv',times_original)
csvwrite('times_biclustering_aa2.csv',times_aa2)
csvwrite('resnorm_biclustering_original.csv',res_norm_hist1)
csvwrite('resnorm_biclustering_aa2.csv',res_norm_hist2)
