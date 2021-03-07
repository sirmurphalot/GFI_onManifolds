function J = Jacobian(mu1,mu2,mu3,dataX,dataY,dataZ, shift_vector)
    % This function calculates the Jacobian value assuming that the
    % constraint is NOT included in the DGA.  
    
    % First, calculate the unconstrained Jacobian matrix
    mat1 = zeros(length(dataX), 3);
    mat1(:,1)=ones(length(dataX),1);
%     mat1(:,4)=dataX-mu1;
    mat2 = zeros(length(dataY), 3);
    mat2(:,2)=ones(length(dataY),1);
%     mat2(:,5)=dataY-mu2;
    mat3 = zeros(length(dataZ), 3);
    mat3(:,3)=ones(length(dataZ),1);
%     mat3(:,6)=dataZ-mu3;
    full_mat=[mat1;mat2;mat3];
    
    % Next, use the constraint that mu1^2+mu2^2+mu3^2=1 to calculate the
    % projection matrix onto the the nullspace of nabla g
    mu1 = mu1 - shift_vector(1);
    mu2 = mu2 - shift_vector(2);
    mu3 = mu3 - shift_vector(3);
    P_mat= eye(3) - ([mu1,mu2,mu3].')*[mu1,mu2,mu3];
    % Perform a singular value decomposition on P_mat
    [U,~,~]=svds(P_mat,2); 
    % Use these columns of U to calculate our projected Jacobian matrix
    J_mat = full_mat*U;
    % Use this to calculate the Jacobian.
    JacValue = -0.5*log(det((J_mat.')*J_mat));
    J = JacValue;
end