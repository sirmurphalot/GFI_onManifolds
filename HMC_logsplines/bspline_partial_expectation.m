function val = bspline_partial_expectation(knots,thetas,spline_number,y_val)
    i = spline_number;
    [d1,d2] = size(thetas);
    d = max(d1,d2);
    if y_val < knots(i)
        val = 0;
    elseif y_val >= knots(i+2)
        val = bspline_expectation(knots,thetas,spline_number);
    elseif (knots(i) < y_val) && (y_val <= knots(i+1))
        if i == 1
            val = 0;
        else
            val = exp(1).^thetas(i-1).*(knots(i)+(-1).*knots(i+1)).*(thetas(i)+(-1).*thetas(i-1)).^(-2).*((-1)+ ...
                    exp(1).^(y_val.*(thetas(i)+(-1).*thetas(i-1))).*(1+y_val.*((-1).*thetas(i)+thetas(i-1))));
        end
    else
        if i == 1
            val = exp(1).^thetas(i).*(knots(i+1)+(-1).*knots(i+2)).*(thetas(i)+(-1).*thetas(i+1)).^(-2).*(1+(-1).*thetas(i)+ ...
                  thetas(i+1)+(-1).*exp(1).^(y_val.*((-1).*thetas(i)+thetas(i+1))).*(1+((-1)+y_val).*thetas(i)+...
                  thetas(i+1)+(-1).*y_val.*thetas(i+1)));
        elseif i == d
            val = bspline_expectation(knots,thetas,spline_number);
        else
            val = (-1).*(knots(i)+(-1).*knots(i+1)).*(exp(1).^thetas(i-1)+exp(1).^thetas(i).*((-1)+thetas(i)+(-1).* ...
                  thetas(i-1))).*(thetas(i)+(-1).*thetas(i-1)).^(-2)+exp(1).^thetas(i).*(knots(i+1)+(-1).*knots(i+2)).*(thetas(i)+( ...
                  -1).*thetas(i+1)).^(-2).*(1+(-1).*thetas(i)+thetas(i+1)+(-1).*exp(1).^(y_val.*((-1).*thetas(i)+ ...
                  thetas(i+1))).*(1+((-1)+y_val).*thetas(i)+thetas(i+1)+(-1).*y_val.*thetas(i+1)));
        end
    end

end