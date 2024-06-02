function bic = compute_BIC(s, s0, lambda, n, d1, d2)
    df = 0;
    q = length(s);
    for i=1:q
        if s(i) == 0
            df = df;
        else
            a = s0(i);
            % compute sum1
            sum1 = 0;
            for j=1:d1
                if j~= i
                    sum1 = sum1 + a*(a-lambda)/(a^2 - s0(j)^2); 
                end
            end
            % compute sum2
            sum2 = 0;
            for j=1:d2
                if j~= i
                    sum2 = sum2 + a*(a-lambda)/(a^2 - s0(j)^2); 
                end
            end
            df = df+1+sum1+sum2;
        end       
    end
    
    bic = n - d1*d2 +log(n)*df;
end
