function obj = obj_SPRR(Y, A, X, lambda1, lambda2)

    s1 = 0.5*norm(Y - A*X, 'fro')^2;
    s2 = lambda1*sum(vecnorm(X'));  % vecnorm(X') returns the 2-norm of each row of X
    s3 = lambda2*sum(svd(X));
    obj = s1+s2+s3;
end
