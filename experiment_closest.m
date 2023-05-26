y = randn(512);
y = (y+y')/2;
tic
[xhat1, res_norm_hist] = dykstra(y, 'closest kinship', 0, 'acceleration', 'original');
toc
tic
[xhat2, res_norm_hist] = dykstra(y, 'closest kinship', 0, 'acceleration', 'aa2');
toc

y = csvread('lung500.csv',1,1);
[xhat, res_norm_hist] = dykstra(y, 'biclustering', 3e5, 'acceleration', 'aa2', 'printevery', 1);