rng('default')
times_original = zeros(10,1);
times_aa2 = zeros(10,1);
for i=1:10
    y = randn(2^i);
    y = (y+y')/2;
    tic
        [xhat1, res_norm_hist1] = dykstra_splitting(y, 'closest kinship', 0, 'acceleration', 'original');
    times_original(i) = toc;
    tic
        [xhat2, res_norm_hist2] = dykstra_splitting(y, 'closest kinship', 0, 'acceleration', 'aa2');
    times_aa2(i) = toc;
end
save("results/closest.mat","res_norm_hist1","res_norm_hist2","times_original","times_aa2")
