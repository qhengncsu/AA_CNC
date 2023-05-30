function [best_lambda] = bootstrap_cvxcluster(U, lambda_seq, B)
n_lambda = length(lambda_seq);
n = size(U,1);
accs = zeros(n_lambda,B);
for i=1:n_lambda
    Xhat = srls_GMC_cvxcluster(U, lambda_seq(i), acceleration='aa2', gamma=0, K=10, phi=1);
    Y = pdist(Xhat);
    Z = linkage(Y);
    T = cluster(Z,"cutoff",0.01,'Criterion','distance');
    if max(T) == 1
        continue
    end
    for j=1:B
        U1 = datasample(U,n);
        U2 = datasample(U,n);
        Xhat1 = srls_GMC_cvxcluster(U1, lambda_seq(i), acceleration='aa2', gamma=0, K=10, phi=1);
        Xhat2 = srls_GMC_cvxcluster(U2, lambda_seq(i), acceleration='aa2', gamma=0, K=10, phi=1);
        Y1 = pdist(Xhat1);
        Z1 = linkage(Y1);
        Y2 = pdist(Xhat2);
        Z2 = linkage(Y2);
        clusters = min(max(T),10);
        T1 = cluster(Z1,'MaxClust',clusters);
        T2 = cluster(Z2,'MaxClust',clusters);
        if max(T1)==1
            RIs(i,j) = 0;
            continue
        else
            U1centers = zeros(max(T1),size(U1,2));
            for k=1:max(T1)
                U1centers(k,:) = mean(U1(T1==k,:));
            end
            Mdl1 = fitcknn(U1centers,1:max(T1),'NumNeighbors',1);
            pred1 = predict(Mdl1,Xhat);
            U2centers = zeros(max(T2),size(U2,2));
            for k=1:max(T2)
                U2centers(k,:) = mean(U2(T2==k,:));
            end
            Mdl2 = fitcknn(U2centers,1:max(T2),'NumNeighbors',1);
            pred2 = predict(Mdl2,Xhat);
            accs(i,j) = cluster_acc(pred1,pred2);
        end
    end
end
end