function J = Jacobian(knots,coef_matrix,thetas,data)
 
    n = length(data);
    d = length(thetas);
    
    % First, calculate the unconstrained Jacobian matrix
    jac_mat = zeros(n,d);
    spline_indices = 1:length(knots);
    num_of_bsplines = length(thetas);
    fun = @(i) (logspline_density(data(i),knots,coef_matrix,thetas));
    density_values = arrayfun(fun,1:n);
    
    fun = @(i) (bspline_expectation(knots,coef_matrix,thetas,i));
    expectation_values = arrayfun(fun,1:d);
    
    
    for i=1:n
        for j=1:d
            knot_interval_number = spline_indices(1,...
                find((data(i) > knots)==1,1,'last'));
            
            if (j>knot_interval_number)
                continue
            elseif (num_of_bsplines+j)<knot_interval_number
                jac_mat(i,j) = -expectation_values(j)/density_values(i);
            else
                spline_indices = 1:length(knots);
                spline_number = spline_indices(1,...
                    find((x_value > knots)==1,1,'last'));
                fun = @(x) (get_bspline_value(x,knots,coef_matrix,spline_number)*...
                    logspline_density(x, knots,coef_matrix,thetas));
                jac_mat(i,j)=integral(fun,knots(1),x_value)/density_values(i);
            end   
        end
    end

    
    % Next, use the constraint that mu1^2+mu2^2+mu3^2=1 to calculate the
    % projection matrix onto the the nullspace of nabla g
    P_mat= eye(d) - expectation_values*(expectation_values')./(sum(expectation_values.^2));
    % Perform a singular value decomposition on P_mat
    [U,~,~]=svds(P_mat,d-1); 
    % Use these columns of U to calculate our projected Jacobian matrix
    J_mat = full_mat*U;
    % Use this to calculate the Jacobian.
    JacValue = 0.5*log(det((J_mat.')*J_mat));
    J = JacValue;
end