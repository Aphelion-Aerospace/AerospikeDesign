function [a,b] = my_bisection(fun,a0,b0,tol)

while (abs(a0-b0) > tol)
    c = (a0+b0)./2;
    if (sign(fun(a0))==sign(fun(c)))
        a0 = c;
    else
        b0 = c;
    end
end
a = a0;
b = b0;