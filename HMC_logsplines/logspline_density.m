function logdens = logspline_density(y_vals,knots,thetas)

    [n,~] = size(y_vals);
    logdens = 0;
    [d1,d2] = size(thetas);
    d = max(d1,d2);
    
    for i=1:n
        for j=1:d
            logdens = logdens + thetas(j)*get_bspline_value(y_vals(i),knots,j);
        end
    end
end