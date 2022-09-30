function ll = ll_density(thetas,knots, data_vals)
    % Function to get the fiducial density for a (x,y) pair.
    % For this implementation, the constraint is NOT built into the DGA.
    % Part of this implementation is projecting the Jacobian onto the
    % constrained space.
    
    % Grab parameters.
    [d1,d2] = size(thetas);
    d = max(d1,d2);
    n = length(data_vals);
    
    expectation_values = zeros(1,d);
    for i=1:d
       expectation_values(i) = bspline_expectation(knots,thetas,i); 
    end
    
    jac_matrix = zeros(n,d);
    for i=1:n
        for j=1:d
            jac_matrix(i,j) = bspline_partial_expectation(knots,thetas,j,data_vals(i))/...
                exp(logspline_density(data_vals(i),knots,thetas));
        end
    end
    
    
    C_vector = expectation_values./sqrt(sum(expectation_values.^2));
    P_mat = diag(d) - C_vector' * C_vector;
    [U,~,~]=svds(P_mat,d-1); 
    % Use these columns of U to calculate our projected Jacobian matrix
    J_mat = jac_matrix*U;
    
    JacValue = 0.5*log(det((J_mat.')*J_mat));
    ll = logspline_density(data_vals,knots,thetas) + JacValue;
end



