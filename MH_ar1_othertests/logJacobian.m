function J = logJacobian(parameters, data_values)
    % This function calculates the Jacobian value assuming that the
    % constraint is NOT included in the DGA.  
    
    % General parameters that will be needed.
    [A,Lambda] = uncollapseParameters(parameters);
    [d,~] = size(Lambda);
    IpA = eye(d) + A;
    IpAinv = inv(IpA);
    
    % Jacobian calculation values from Murph et al
    %%% I'll assume that cols are timepoints
    n_time = length(data_values);
    
    %%% The following excludes mu -- I'm setting it equal to zero.
    J_mat = zeros(n_time, (d-1)*d/2 + d);
    count = 0;
    lower_index = 1;
    upper_index = n_time;
    
    % Calculate A partials
    % The A partials are calculated BY ROW.
    for i=1:(d-1)
        for j=(i+1):d
            Jij = zeros(d,d);
            Jij(i,j) = 1;
            count = count + 1;
            J_mat(1:n_time,count)=2*(IpAinv)*...
                (-Jij+Jij')*(IpAinv')*(data_values');
        end
    end
    % Calculate Lambda partials
    for i=1:d
        Jii = zeros(d,d);
        Jii(i,i) = 1;
        count = count + 1;
        J_mat(1:n_time,count)=(Lambda(i,i).^(-1))*(IpA')*...
            (IpAinv)*Jii*((IpAinv')*IpA)*(data_values');
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
    JacValue = 0.5*sum(log(eig((J_mat.')*J_mat)));
    J = JacValue;
end