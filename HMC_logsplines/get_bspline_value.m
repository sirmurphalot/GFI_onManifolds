function val = get_bspline_value(x_value,knots,coef_matrix,bspline_number)
    [~,bspline_dim]=size(coef_matrix);
    indices = bspline_number:(bspline_dim+bspline_number);
    spline_intervals = knots(indices);
    
    bsplines_coefs=coef_matrix(((bspline_number-...
        1)*bspline_dim+1):(bspline_number)*bspline_dim,:);
    
    powers = 0:(bspline_dim-1);
    x_poly_terms = flip(arrayfun(@(i)(x_value^i),powers))';
    
    val = 0;
    for i=2:length(spline_intervals)
        if (x_value >= spline_intervals(i-1))&&(x_value < spline_intervals(i))
            val = bsplines_coefs(i-1,:)*x_poly_terms;
        end
    end
end