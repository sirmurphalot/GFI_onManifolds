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
    [n_obs,n_time] = size(data_values);
    
    %%% The following excludes mu -- I'm setting it equal to zero.
    J_mat = zeros(n_time*n_obs, (d-1)*d/2 + d);
    count = 0;
    lower_index = 1;
    upper_index = n_time;
    
    % Calculate A partials
    % The A partials are calculated BY ROW.
    obs_count = 0;
    for n = 1:n_obs
        curr_obs = data_values(n,:);
        count = 0;
        for i=1:(d-1)
            for j=(i+1):d
                Jij = zeros(d,d);
                    Jij(i,j) = 1;
                    count = count + 1;
                zz = (1:n_time)+(obs_count*n_time);
                J_mat(zz,count)=2*(IpAinv)*...
                        (-Jij+Jij')*(IpAinv')*(curr_obs');
            end
        end

        % Calculate Lambda partials
        for i=1:d
            Jii = zeros(d,d);
            Jii(i,i) = 1;
            count = count + 1;
            zz = (1:n_time)+(obs_count*n_time);
            J_mat(zz,count)=(Lambda(i,i).^(-1))*(IpA')*...
                            (IpAinv)*Jii*(IpA*(IpAinv'))*(curr_obs');
        end
        obs_count = obs_count + 1;
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
