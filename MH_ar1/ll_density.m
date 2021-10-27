function ll = ll_density(parameters, data_values)
    [n,d] = size(data_values);
    logDens=0;
    [A,Lambda] = uncollapseParameters(parameters);
    U_mat = (eye(d)-A)/(eye(d)+A);
    Lambda_inv = diag(diag(Lambda).^(-2));
    precisionMatrix = U_mat * Lambda_inv * U_mat';
    
    for i=1:n
        % Typical normal density assuming data is iid (w/in groups,
        % just independent between groups).
        logDens=logDens-(d/2)*log(2*pi)+(1/2)*log(det(precisionMatrix))-...
            (1/2)*(data_values(i,:)')'*precisionMatrix*(data_values(i,:)'); 
    end
    
    
    jac = Jacobian(parameters, data_values);
    
   ll = logDens+jac;
end



