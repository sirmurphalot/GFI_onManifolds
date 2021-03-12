function dens = logspline_density(x_value,knots,coef_matrix,thetas)
    spline_indices = 1:length(knots);
    knot_interval_number = spline_indices(1,...
        find((x_value > knots)==1,1,'last'));
    [~,bspline_dim] = size(coef_matrix);
    num_of_bsplines = length(thetas);
    total_num_of_knots = length(knots);
    indices = 1:num_of_bsplines;
    
    which_work = ((indices + bspline_dim)>=knot_interval_number)&...
        ((indices)<=total_num_of_knots);
    fun = @(i) (get_bspline_value(x_value,knots,...
        coef_matrix,i)*thetas(i));
    dens = exp(sum(arrayfun(@(i)(fun(i)),which_work)));
end