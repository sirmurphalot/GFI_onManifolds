function [J,dJ] = Jacobian(knots,thetas,data)
 
    n = length(data);
    d = length(thetas);
    
    % First, calculate the unconstrained Jacobian matrix
    jac_mat = zeros(n,d);
    
    % To cut down on excess computations, we'll precalculate the densities:
    fun = @(i) (exp(logspline_density(data(i),knots,thetas)));
    density_values = arrayfun(fun,1:n);
    % To cut down on excess computations, we'll precalculate the expectations:
    fun = @(i) (bspline_expectation(knots,thetas,i));
    expectation_values = arrayfun(fun,1:d);
    expectation_norm = sum(expectation_values.^2);
    % To cut down on excess computations, we'll precalculate the 2nd degree
    % expectations:
    expectations2nd_values = zeros(d,d);
    for i=1:d
        for j=1:d
           expectations2nd_values(i,j)=bspline_2nd_expectation(knots,...
               thetas,i,j);
        end
    end
    
    % Two third-order arrays we'll need to fill in for the derivative
    % values:
    d_jac_mat = zeros(n,d,d);
    d_P_mat = zeros(d,d,d);
    dC_mat = zeros(1,d,d);
    
    % Fill in the jacobian matrix and derivative of jacobian matrix
    for i=1:n
        for j=1:d

            jac_mat(i,j) = -bspline_partial_expectation(knots,thetas,j,data(i))/density_values(i);
            for q=1:d
                d_jac_mat(i,j,q) = (bspline_2nd_partial_expectation(knots,thetas,...
                    j,q,data(i)) - get_bspline_value(data(i),knots,q)*...
                    bspline_partial_expectation(knots,thetas,j,data(i)))/density_values(i);
            end  
        end
    end

    for i=1:d
        for j=1:d
            dC_mat(1,j,i)=(expectation_norm)^(-3/2) * (-expectation_values(i)*...
                (expectation_values * expectations2nd_values(:,j)) + ...
                expectation_norm * expectations2nd_values(i,j));
        end
    end
    
    % Now, calculate the regular projection matrix:
    full_mat = jac_mat;
    P_mat= eye(d) - (expectation_values'*expectation_values)./(expectation_norm);
    % Perform a singular value decomposition on P_mat
    try
        [U,~,~]=svds(P_mat,d-1); 
    catch
        disp("hey");
    end
    % Use these columns of U to calculate our projected Jacobian matrix
    J_mat = full_mat*U;
    % Use this to calculate the Jacobian.
    JacValue = 0.5*log(det((J_mat.')*J_mat));
    J = JacValue;
    
    dJ=zeros(1,d);
    expectation_values = expectation_values./sqrt((expectation_norm));
    for q=1:d
        d_P_mat=-(dC_mat(1,:,q)')*(expectation_values)-...
            (expectation_values')*dC_mat(1,:,q);
        d_J_mat = d_jac_mat(:,:,q);
        penrose_inv = pinv(full_mat*P_mat*full_mat');
        dJ(q)=0.5*trace(penrose_inv*(d_J_mat*P_mat*full_mat'+...
            full_mat*d_P_mat*full_mat'+full_mat*P_mat*d_J_mat'));
    end
    
end