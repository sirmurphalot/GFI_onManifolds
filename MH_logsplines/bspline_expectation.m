function val = bspline_expectation(knots,thetas,spline_number)
    i = spline_number;
    [d1,d2] = size(thetas);
    d = max(d1,d2);
    if i == 1
        val = (knots(i+1)-knots(i))*(1 + exp(thetas(i)) * (thetas(i)-1)) / thetas(i)^2;
        val = val + (knots(i+1) - knots(i+2)) * (-exp(thetas(i+1)) + exp(thetas(i)) * (thetas(i) - 1 - thetas(i+1) )) / (thetas(i) - thetas(i+1))^2;
    elseif i == d
        val = (knots(i+1) - knots(i)) * (-exp(thetas(i-1)) + exp(thetas(i)) * (thetas(i)-1-thetas(i-1))) / (thetas(i)-thetas(i-1))^2;
        val = val + (knots(i+2) - knots(i+1)) * (1 + exp(thetas(i)) * (thetas(i)-1)) / thetas(i)^2;
    else
        val = (knots(i+1) - knots(i)) * (-exp(thetas(i-1)) + exp(thetas(i)) * (thetas(i)-1-thetas(i-1))) / (thetas(i)-thetas(i-1))^2;
        val = val + (knots(i+1) - knots(i+2)) * (-exp(thetas(i+1)) + exp(thetas(i)) * (thetas(i) - 1 - thetas(i+1) )) / (thetas(i) - thetas(i+1))^2;
    end
end