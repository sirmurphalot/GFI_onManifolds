function [c, ceq] = optim_constraint(knots, thetas)
    ceq = logspline_constraint(knots, thetas);
    c = [];
end