function val = bspline_2nd_expectation(knots,thetas,spline_number1,...
    spline_number2)

    splinenumber1 = min(spline_number1, spline_number2);
    splinenumber2 = max(spline_number1, spline_number2); 
    i = splinenumber1;
    [d1,d2] = size(thetas);
    d = max(d1,d2);
    
    if (abs(splinenumber1 - splinenumber2) > 1)
        val = 0;
    elseif (abs(splinenumber1 - splinenumber2) == 1)
        val = (knots(i+1)+(-1).*knots(i+2)).*(thetas(i)+(-1).*thetas(i+1)).^(-3).*(exp(1).^thetas(i+1).*((-2)+...
                (-1).*thetas(i)+thetas(i+1))+exp(1).^thetas(i).*(2+(-1).*thetas(i)+thetas(i+1)));
    else
        if splinenumber1 == 1
            val = (-1).*(knots(i+1)+(-1).*knots(i+2)).*(thetas(i)+(-1).*thetas(i+1)).^(-3).*((-2).*...
                exp(1).^thetas(i+1)+exp(1).^thetas(i).*(2+thetas(i).^2+2.*thetas(i+1)+...
                thetas(i+1).^2+(-2).*thetas(i).*(1+thetas(i+1))));
        elseif splinenumber2 == d
            val = (-1).*(knots(i)+(-1).*knots(i+1)).*(thetas(i)+(-1).*thetas(i-1)).^(-3).*((-2).*exp(1).^thetas(i-1)+ ...
                    exp(1).^thetas(i).*(2+thetas(i).^2+2.*thetas(i-1)+thetas(i-1).^2+(-2).*thetas(i).*(1+thetas(i-1))));
        else
            val = (-1).*(knots(i)+(-1).*knots(i+1)).*(thetas(i)+(-1).*thetas(i-1)).^(-3).*((-2).*exp(1).^thetas(i-1)+ ...
                  exp(1).^thetas(i).*(2+thetas(i).^2+2.*thetas(i-1)+thetas(i-1).^2+(-2).*thetas(i).*(1+thetas(i-1))))+(-1).*( ...
                  knots(i+1)+(-1).*knots(i+2)).*(thetas(i)+(-1).*thetas(i+1)).^(-3).*((-2).*exp(1).^thetas(i+1)+exp( ...
                  1).^thetas(i).*(2+thetas(i).^2+2.*thetas(i+1)+thetas(i+1).^2+(-2).*thetas(i).*(1+thetas(i+1))));
        end      
    end
end