function dens = logspline_density(x_value,knots,coef_matrix,thetas)
    spline_indices = 1:length(knots);
    [~,bspline_dim] = size(coef_matrix);
    num_of_bsplines = length(thetas);
    total_num_of_knots = length(knots);
    indices = 1:num_of_bsplines;
    
    [n,d] = size(x_value);
    dens = zeros(size(n,d));
    
    for i=1:n
        for j=1:d
            knot_interval_number = spline_indices(1,...
                find((x_value(i,j) > knots)==1,1,'last'));
            which_work = ((indices + bspline_dim)>=knot_interval_number)&...
                ((indices)<=total_num_of_knots);
            fun = @(k) (get_bspline_value(x_value(i,j),knots,...
                coef_matrix,k)*thetas(k));
            dens(i,j) = exp(sum(arrayfun(@(k)(fun(k)),indices(which_work))));
        end
    end
end