function y = doFilter(x,a)
    N = length(x);
    b = 1 - a;

    y = zeros(1, N);

    y(1) = x(1);

    for i = 2:N
        y(i) = a*y(i-1) + b*x(i);
    end
end

