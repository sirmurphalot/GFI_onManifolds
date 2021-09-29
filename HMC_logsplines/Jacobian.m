function [J,dJ] = Jacobian(knots,coef_matrix,thetas,data)
 
    n = length(data);
    d = length(thetas);
    [~,bspline_dim] = size(coef_matrix);
    
    % First, calculate the unconstrained Jacobian matrix
    jac_mat = zeros(n,d);
    spline_indices = 1:(length(knots)-2);
    num_of_bsplines = length(thetas);
    
    % Two third-order arrays we'll need to fill in for the derivative
    % values:
    d_jac_mat = zeros(n,d,d);
    d_P_mat = zeros(d,d,d);
    dC_mat = zeros(1,d,d);
    
    % Fill in the jacobian matrix and derivative of jacobian matrix
    for i=1:n
        for j=1:d
            knot_interval_number = spline_indices(1,...
                find((data(i) > knots)==1,1,'last'));
            
            if (j>knot_interval_number)
                continue
            elseif (num_of_bsplines+j)<knot_interval_number
                jac_mat(i,j) = -expectation_values(j)/density_values(i);
                for q=1:d
                    if (q>knot_interval_number)
                        continue
                    elseif (num_of_bsplines+q)<knot_interval_number
                        d_jac_mat(i,j,q) = ...
                            -(expectations2nd_values(j,q))/density_values(i)-...
                            jac_mat(i,j)*...
                            get_bspline_value(data(i),knots,coef_matrix,q);
                    end
                end
            else
                fun = @(x) (get_bspline_value(x,knots,coef_matrix,j)*...
                    logspline_density(x, knots,coef_matrix,thetas));
                jac_mat(i,j)=integral(fun,knots(1),data(i))/density_values(i);
                for q=1:d
                    if (q>knot_interval_number)
                        continue
                    else
                        fun = @(x) (get_bspline_value(x,knots,coef_matrix,j).*...
                            logspline_density(x,knots,coef_matrix,thetas).*...
                            get_bspline_value(x,knots,coef_matrix,q));
                        val_2ndegree=integral(fun,knots(1),data(i))/density_values(i);
                        d_jac_mat(i,j,q)=-val_2ndegree + ...
                            jac_mat(i,j)*...
                            get_bspline_value(data(i),knots,coef_matrix,q);
                    end
                end
                
            end   
        end
    end

    for i=1:d
        for j=1:d
            dC_mat(1,j,i)=(expectation_norm*expectations2nd_values(i,j)-...
                expectation_values(i)*(sum(expectations2nd_values(i,:))))/...
                (expectation_norm^(3/2));
        end
    end
    
    % Now, calculate the regular projection matrix:
    full_mat = jac_mat;
    P_mat= eye(d) - (expectation_values'*expectation_values)./(expectation_norm);
    % Perform a singular value decomposition on P_mat
    [U,~,~]=svds(P_mat,d-1); 
    % Use these columns of U to calculate our projected Jacobian matrix
    J_mat = full_mat*U;
    % Use this to calculate the Jacobian.
    JacValue = 0.5*log(det((J_mat.')*J_mat));
    J = JacValue;
    
    dJ=zeros(1,d);
    for q=1:d
        d_P_mat=-(dC_mat(1,:,q)')*(expectation_values)-...
            (expectation_values')*dC_mat(1,:,q);
        d_J_mat = d_jac_mat(:,:,q);
        penrose_inv = pinv(full_mat*P_mat*full_mat');
        dJ(q)=0.5*trace(penrose_inv*(d_J_mat*P_mat*full_mat'+...
            full_mat*d_P_mat*full_mat'+full_mat*P_mat*d_J_mat'));
    end
    
end