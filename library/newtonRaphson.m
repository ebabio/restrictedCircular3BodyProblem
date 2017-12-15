function x = newtonRaphson(f, dfdx, x0, e)

if(nargin ==3)
    e = 1e-12;
else
    e = abs(e);
end

error = 2*e;
x = x0;
counter = 0;
while(abs(error)>e)
    error = f(x)/dfdx(x);
    x = x - error;
    counter = counter +1;
    if(counter > 3e1)
        x = NaN;
        break
    end
end
