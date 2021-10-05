function c = logspline_constraint(knots, thetas)
    [~, knot_length] = size(knots);
    c = 0;
    c = c + (exp(thetas(1)) - 1)*((knots(2) - knots(1))/(thetas(1)));
    for i=2:(knot_length-2)
        c = c + (exp(thetas(i)) - exp(thetas(i-1)))*((knots(i+1) - knots(i))/(thetas(i) - thetas(i-1)));
    end
    c = c + (exp(thetas(knot_length-2)) - 1)*((knots(knot_length) - knots(knot_length-1))/(thetas(knot_length-2)));
    c = log(c);
end