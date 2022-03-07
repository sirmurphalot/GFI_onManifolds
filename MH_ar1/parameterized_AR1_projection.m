function [distance_from_proposal] = parameterized_AR1_projection(rho,...
    sd2, Tx, curr_loc_x, proposal_value, n_timepoints)
    
    % Start by getting the covariance matrix from the given rho and sigma.
    covariance_matrix = zeros(n_timepoints, n_timepoints);
    for i=1:n_timepoints
        for j=1:n_timepoints
            covariance_matrix(i,j) = rho^(abs(i-j));
        end
    end
    covariance_matrix = covariance_matrix.*((sd2)/(1-rho.^2));
    [U, Lambda2] = eig(covariance_matrix);
    Lambda = diag(diag(Lambda2.^(0.5)));
    
    % We will attempt to match signature matrices of our eigendecomposition 
    % with the eigendecomposition we did to get Tx.
    
    % There are three values.  The ones from the sd2,rho; the ones from the
    % current location; the ones from the proposal value.  We will gather
    % the Lambdas and Us from each and ALIGN their signature matrices.
    
    %%% This is aligning the values from the sd2,rho.
    [sorted_Lambda, sorted_indices] = sort(diag(Lambda));
    U_sorted = U(:,sorted_indices);
    
    %%% This is aligning the values for the current location.
    [A_x, Lambda_x] = uncollapseParameters(curr_loc_x);
    [sorted_Lambda_x, sorted_indices_x] = sort(diag(Lambda_x));
    Ux = (eye(n_timepoints) - A_x)/(eye(n_timepoints) + A_x);
    Ux_sorted = Ux(:,sorted_indices_x);
    Ax_sorted = (eye(n_timepoints)-Ux_sorted)/(eye(n_timepoints)+Ux_sorted);
    curr_location_sorted = collapseParameters(Ax_sorted, diag(sorted_Lambda_x));
    
    %%% This is aligning the values for the proposal.
    [A_pro, Lambda_pro] = uncollapseParameters(proposal_value);
    [sorted_Lambda_pro, sorted_indices_pro] = sort(diag(Lambda_pro));
    U_pro = (eye(n_timepoints) - A_pro)/(eye(n_timepoints) + A_pro);
    U_pro_sorted = U_pro(:,sorted_indices_pro);
    A_pro_sorted = (eye(n_timepoints)-U_pro_sorted)/(eye(n_timepoints)+U_pro_sorted);
    proposal_value_sorted = collapseParameters(A_pro_sorted, diag(sorted_Lambda_pro));
    
    %%% Now we determine a signature matrix for the sd2,rho value based on
    %%% its correlation with the eigenvecters for the current location.
    [~,n_vectors] = size(Ux_sorted);
    signature_matrix = zeros(n_vectors,1);
    for k = 1:n_vectors
        temp_sign = sign(corr(U_sorted(:,k), Ux_sorted(:,k)));
        signature_matrix(k) = 1*temp_sign;
    end
    signature_matrix = diag(signature_matrix);
    
    U_final = U_sorted * signature_matrix;
%     E = eig(U_final);
    
%     if ~(sum(abs(real(E)+1)<1e-4)==0)
%         distance_from_proposal = 1e4;
%     else
    I = eye(n_timepoints);
    A_final = (I-U_final)/(I+U_final);

    % We will shift this parameter value associated with rho,sd2 by the
    % current location, project it onto the tangent plance, and shift
    % back.
    q_shifted = collapseParameters(A_final, diag(sorted_Lambda)) - ...
        curr_location_sorted;
    projection_matrix = (Tx/((Tx')*Tx))*(Tx');

    value_on_tangent_space_nonshifted = projection_matrix*q_shifted;
    value_on_tangent_space = value_on_tangent_space_nonshifted + ...
        curr_location_sorted;

    % Finally, this method returns the distance between the sd2,rho
    % value and the proposal value.
    distance_from_proposal = sqrt((value_on_tangent_space')*...
        proposal_value_sorted);
%     end
end