function [c, dc] = logspline_constraint(knots, thetas)
    [knot_length,~] = size(knots);
    [d,~] = size(thetas);
    c = 0;
    c = c + (exp(thetas(1)) - 1)*((knots(i) - knots(i+1))/(theta(i)));
    for i=2:(knot_length-2)
        c = c + (exp(thetas(i)) - exp(thetas(i-1)))*((knots(i+1) - knots(i))/(theta(i) - theta(i-1)));
    end
    c = c + (exp(thetas(knot_length-2)) - 1)*((knots(knot_length) - knots(knot_length-1))/(theta(knot_length-2)));
    c = log(c);
    
    % Assuming c=0, we can cut back a bit on computations.
    dc = zeros(d);
    for i=1:d
        dc(i) = bspline_expectation(knots,thetas,i);
    end
end