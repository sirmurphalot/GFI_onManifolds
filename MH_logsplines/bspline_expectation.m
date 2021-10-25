function val = bspline_expectation(knots,thetas,spline_number)
    i = spline_number;
    [d1,d2] = size(thetas);
    d = max(d1,d2);
    if i == 1
        val = (knots(i+1)+(-1).*knots(i+2)).*(thetas(i)+(-1).*thetas(i+1)).^(-2).*((-1).*exp(1).^thetas(i+1)+exp( ...
                1).^thetas(i).*(1+(-1).*thetas(i)+thetas(i+1)));
    elseif i == d
        val = (-1).*(knots(i)+(-1).*knots(i+1)).*(exp(1).^thetas(i-1)+exp(1).^thetas(i).*...
                ((-1)+thetas(i)+(-1).* ...
                thetas(i-1))).*(thetas(i)+(-1).*thetas(i-1)).^(-2);
    else
        val = (-1).*(knots(i)+(-1).*knots(i+1)).*(exp(1).^thetas(i-1)+exp(1).^thetas(i).*((-1)+thetas(i)+(-1).* ...
                  thetas(i-1))).*(thetas(i)+(-1).*thetas(i-1)).^(-2)+(knots(i+1)+(-1).*knots(i+2)).*(thetas(i)+(-1).*...
                  thetas(i+1)).^(-2).*((-1).*exp(1).^thetas(i+1)+exp(1).^thetas(i).*(1+(-1).*thetas(i)+thetas(i+1)));
    end
    if isinf(val)
        disp("hey");
    end
end