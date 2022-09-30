function val = bspline_2nd_partial_expectation(knots,thetas,spline_number1,...
    spline_number2, y_val)

    splinenumber1 = min(spline_number1, spline_number2);
    splinenumber2 = max(spline_number1, spline_number2); 
    i = splinenumber1;
    [d1,d2] = size(thetas);
    d = max(d1,d2);
    
    if (abs(splinenumber1 - splinenumber2) > 1)
        val = 0;
    elseif (abs(splinenumber1 - splinenumber2) == 1)
        if y_val >= knots(splinenumber2+1)
            val = bspline_2nd_expectation(knots,thetas,spline_number1,spline_number2);
        elseif y_val < knots(splinenumber2)
            val = 0;
        else
            val = exp(1).^thetas(i).*(knots(i+1)+(-1).*knots(i+2)).*(thetas(i)+(-1).*thetas(i+1)).^(-3).*(2+(-1).*thetas(i)+ ...
                  thetas(i+1)+(-1).*exp(1).^(y_val.*((-1).*thetas(i)+thetas(i+1))).*(2+((-1)+y_val).*y_val.*thetas(i).^2+ ...
                  thetas(i+1)+y_val.*thetas(i+1).*((-2)+((-1)+y_val).*thetas(i+1))+thetas(i).*((-1)+2.*y_val.*(1+thetas(i+1)+(-1).* ...
                  y_val.*thetas(i+1)))));
        end
    else
        if y_val >= knots(splinenumber2+2)
            val = bspline_2nd_expectation(knots,thetas,spline_number1,spline_number2);
        elseif y_val < knots(splinenumber1)
            val = 0;
        elseif (knots(splinenumber1) <= y_val) && (y_val < knots(splinenumber1+1))
            if splinenumber1 == 1
                val = 0;
            else
                val=exp(1).^thetas(i-1).*(knots(i)+(-1).*knots(i+1)).*((-2)+exp(1).^(y_val.*(thetas(i)+(-1).*thetas(i-1))) ...
                      .*(2+y_val.*((-2)+y_val.*(thetas(i)+(-1).*thetas(i-1))).*(thetas(i)+(-1).*thetas(i-1)))).*((-1).*thetas(i)+ ...
                      thetas(i-1)).^(-3);
            end
        else
            if splinenumber1 == 1
                val = (-1).*exp(1).^thetas(i).*(knots(i+1)+(-1).*knots(i+2)).*(thetas(i)+(-1).*thetas(i+1)).^(-3).*(2+ ...
                      thetas(i).^2+(-2).*thetas(i).*(1+thetas(i+1))+thetas(i+1).*(2+thetas(i+1))+(-1).*exp(1).^(y_val.*((-1).* ...
                      thetas(i)+thetas(i+1))).*(2+((-1)+y_val).^2.*thetas(i).^2+((-1)+y_val).*thetas(i+1).*((-2)+((-1)+y_val).* ...
                      thetas(i+1))+(-2).*((-1)+y_val).*thetas(i).*((-1)+((-1)+y_val).*thetas(i+1))));
            elseif splinenumber2 == d
                val = bspline_2nd_expectation(knots,thetas,spline_number1,spline_number2);
            else
                val = (-1).*(knots(i)+(-1).*knots(i+1)).*(thetas(i)+(-1).*thetas(i-1)).^(-3).*((-2).*exp(1).^thetas(i-1)+ ...
                      exp(1).^thetas(i).*(2+thetas(i).^2+2.*thetas(i-1)+thetas(i-1).^2+(-2).*thetas(i).*(1+thetas(i-1))))+(-1).* ...
                      exp(1).^thetas(i).*(knots(i+1)+(-1).*knots(i+2)).*(thetas(i)+(-1).*thetas(i+1)).^(-3).*(2+thetas(i).^2+( ...
                      -2).*thetas(i).*(1+thetas(i+1))+thetas(i+1).*(2+thetas(i+1))+(-1).*exp(1).^(y_val.*((-1).*thetas(i)+thetas(i+1)) ...
                      ).*(2+((-1)+y_val).^2.*thetas(i).^2+((-1)+y_val).*thetas(i+1).*((-2)+((-1)+y_val).*thetas(i+1))+( ...
                      -2).*((-1)+y_val).*thetas(i).*((-1)+((-1)+y_val).*thetas(i+1))));
            end
        end
    end
end