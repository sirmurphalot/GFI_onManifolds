function logdens = logspline_density(y_vals,knots,thetas)

    [n,~] = size(x_value);
    logdens = 0;
    [d,~] = size(thetas);
    
    for i=1:n
        for j=1:d
            dens = dens + get_bspline_value(y_vals(i),knots,j);
        end
    end
end