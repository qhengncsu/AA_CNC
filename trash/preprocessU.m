function [U, D, w] = preprocessU(U,K,phi)
    if size(U,2)>10
        Idx = knnsearch(U,U,Distance="fasteuclidean",K=K+1);
    else
        Idx = knnsearch(U,U,Distance="euclidean",K=K+1);
    end
    Idx = Idx(:,2:end);
    l = 1;
    pairs = [0,0];
    for i=1:size(Idx,1)
        for j=1:size(Idx,2)
            if i<Idx(i,j)
                pairs(l,:) = [i,Idx(i,j)];
                l = l+1;
            else
                if ~ismember(i,Idx(Idx(i,j),:))
                    pairs(l,:) = [Idx(i,j),i];
                    l = l+1;
                end
            end   
        end
    end
    num_pairs = size(pairs,1);
    D1 = sparse(1:num_pairs,pairs(:,1),ones(num_pairs,1),num_pairs,size(U,1));
    D2 = sparse(1:num_pairs,pairs(:,2),ones(num_pairs,1),num_pairs,size(U,1));
    D = D1-D2;
    w = zeros(num_pairs,1);
    for l=1:num_pairs
        w(l) = exp(-phi*norm(U(pairs(l,1),:)-U(pairs(l,2),:))^2);
    end
end