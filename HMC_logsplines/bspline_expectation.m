function val = bspline_expectation(knots,coef_matrix,thetas,spline_number)
    [~,bspline_dim] = size(coef_matrix);
    fun = @(x)(get_bspline_value(x,knots,coef_matrix,spline_number)*...
        logspline_density(x, knots,coef_matrix,thetas));
    val = integral(fun, knots(spline_number), ...
        knots(spline_number+bspline_dim));
end