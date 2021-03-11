function dens = logspline_density(x_value,knots,coef_matrix,thetas)
    spline_indices = 1:length(knots);
    spline_number = spline_indices(1,find((x_value > knots)==1,1,'last'));
    
    [~,powers] = size(coef_matrix);
    powers=powers-1;
    dens = exp(sum(arrayfun(@(i)(x^i),powers).*...
        coef_matrix(spline_number,:))*...
        thetas(spline_number));
end