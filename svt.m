function x_svt = svt(x,lambda)
    [U,S,V] = svd(x);
    S = max(S-lambda,0);
    x_svt = U*S*V';
end