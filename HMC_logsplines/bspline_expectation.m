function val = bspline_expectation(knots,coef_matrix,thetas,spline_number)
    [~,powers] = size(coef_matrix);
    powers=powers-1;
    
    fun = @(x)(sum(arrayfun(@(i)(x^i),powers).*...
        coef_matrix(spline_number,:))*...
        logspline_density(x, knots,coef_matrix,thetas));
    val = integral(fun, knots(spline_number), knots(spline_number)+1);
end