function J_par = JacobianParameterized(theta,phi,dataX,dataY,dataZ)
    % This function calculates the Jacobian value assuming that the
    % constraint IS included in the DGA. 
    
    % Calculate the Jacobian matrix as one always does for the fiducial
    % approach
    mat1 = zeros(length(dataX), 5);
    mat1(:,1)=(cos(theta)*cos(phi)).*ones(length(dataX),1);
    mat1(:,2)=(-sin(theta)*sin(phi)).*ones(length(dataX),1);
    mat1(:,3)=dataX-(cos(theta)*sin(phi)).*ones(length(dataX),1);
    mat2 = zeros(length(dataY), 5);
    mat2(:,1)=(sin(theta)*cos(phi)).*ones(length(dataX),1);
    mat2(:,2)=(cos(theta)*sin(phi)).*ones(length(dataX),1);
    mat2(:,4)=dataY-(sin(theta)*sin(phi)).*ones(length(dataX),1);
    mat3 = zeros(length(dataZ), 5);
    mat3(:,1)=(-sin(phi)).*ones(length(dataX),1);
    mat3(:,5)=dataZ-cos(phi)*ones(length(dataX),1)-ones(length(dataX),1);
    full_mat=[mat1;mat2;mat3];
    
    % Return the Jacobian value
    J_par=-0.5*log(det((full_mat.')*full_mat));
end
