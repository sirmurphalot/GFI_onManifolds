function J = Jacobian(parameters, data_values)
    % This function calculates the Jacobian value assuming that the
    % constraint is NOT included in the DGA.  
    
    % General parameters that will be needed.
    [A,Lambda] = uncollapseParameters(parameters);
    [d,~] = size(Lambda);
    IpA = eye(d) + A;
    IpAinv = inv(IpA);
    
    % Jacobian calculation values from Murph et al
    %%% I'll assume that rows are observations and cols are timepoints
    [n_obs, n_time] = size(data_values);
    J_mat = zeros(n_obs*n_time, (d-1)*d/2 + d);
    for row=1:n_obs
        count = 0;
        lower_index = ((row-1)*n_time+1);
        upper_index = row*n_time;
        % Calculate A partials
        for i=1:(d-1)
            for j=(i+1):d
                Jij = zeros(d,d);
                Jij(i,j) = 1;
                count = count + 1;
                J_mat(lower_index:upper_index,count)=2*IpAinv*...
                    (-Jij+Jij)*IpAinv*(data_values(row,:)');
            end
        end
        % Calculate Lambda partials
        for i=1:d
            Jii = zeros(d,d);
            Jii(i,i) = 1;
            count = count + 1;
            J_mat(lower_index:upper_index,count)=(Lambda(i,i).^(-1))*(IpA')*...
                IpAinv*Jii*(((IpA')*IpAinv)')*(data_values(row,:)');
        end
    end
    
    % Calculation of P matrix derivative values
    dc = dar1_constraint(parameters);
    
    % Calculate the P matrix, and decompose.
    [~,d2] = size(dc);
    P_mat= eye(d2) - ((dc')/(dc*(dc')))*dc;
    % Perform a singular value decomposition on P_mat
    [U,~,~]=svds(P_mat,2); 
    % Use these columns of U to calculate our projected Jacobian matrix
    J_mat = J_mat*U;
    % Use this to calculate the Jacobian.
    JacValue = -0.5*log(det((J_mat.')*J_mat));
    J = JacValue;
end