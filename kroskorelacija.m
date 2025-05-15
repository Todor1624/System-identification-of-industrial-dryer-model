function rxy = kroskorelacija(x,y,tau)
    
    rxy = zeros(1,2*tau+1);
    N = length(x);

    for i = 0:tau
        ind_pos = tau + 1 + i;
        ind_neg = tau + 1 - i;
        for k = 1 : N-i
            rxy(ind_pos) = rxy(ind_pos) + y(k)*x(k+i);
        end
            rxy(ind_pos) = rxy(ind_pos)/N;

         for k = 1 : N-i
            rxy(ind_neg) = rxy(ind_neg) + x(k)*y(k+i);
         end
            rxy(ind_neg) = rxy(ind_neg)/N;



    end
end

