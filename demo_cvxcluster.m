rng('default') % For reproducibility
U = csvread("movement_libras.csv");
labels = U(:,end);
U = U(:,1:(end-1));
labels_to_filter = [1:15];
rows_to_filter = ismember(labels, labels_to_filter);
U = U(rows_to_filter,:);
labels = labels(rows_to_filter);
unique_labels = unique(labels);
n_labels = length(unique_labels);
label_map = containers.Map(unique_labels, 1:n_labels);
labels = arrayfun(@(x) label_map(x), labels);
ARs_convex = zeros(20,1);
i = 1
for lambda=1:1:20
    Xhat = srls_GMC_cvxcluster(U, lambda, acceleration='aa2', gamma=0, K=5, phi=0.5);
    Y = pdist(Xhat);
    Z = linkage(Y);
    %dendrogram(Z)
    T = cluster(Z,"cutoff",0.01,'Criterion','distance');
    ARs_convex(i) = RandIndex(T,labels);
    i = i+1;
end

ARs = zeros(20,1);
i = 1
for lambda=1:1:20
    Xhat = srls_GMC_cvxcluster(U, lambda, acceleration='aa2', gamma=0.8, K=5, phi=0.5);
    Y = pdist(Xhat);
    Z = linkage(Y);
    %dendrogram(Z)
    T = cluster(Z,"cutoff",0.01,'Criterion','distance');
    ARs(i) = RandIndex(T,labels);
    i = i+1;
end

[Xhat, Vhat, res_norm_hist] = srls_GMC_cvxcluster(U, 50000, acceleration='aa2', gamma=0, K=10, phi=0.5);
scatter(Xhat(:,1),Xhat(:,2))
Y = pdist(Xhat);
Z = linkage(Y);
T = cluster(Z,"cutoff",0.3,'Criterion','distance');
AR = RandIndex(T,label)