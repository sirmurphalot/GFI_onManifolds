function [c, dc] = logspline_constraint(knots,coef_matrix,thetas)
    fun = @(x) (logspline_density(x,knots,coef_matrix,thetas));
    marginal = integrate(fun, 0, 1);
    c = log(marginal);
    
    % Assuming c=0, we can cut back a bit on computations.
    fun = @(i) (bspline_expectation(knots,coef_matrix,thetas,i));
    dc = arrayfun(fun, 1:length(knots))';
end