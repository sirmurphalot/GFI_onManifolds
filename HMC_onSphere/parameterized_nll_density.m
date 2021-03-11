function [nll,dnll] = parameterized_nll_density(q,s1,s2,s3,dataX,dataY,dataZ)
    % Function to get the fiducial density for a (x,y) pair.
    % For this implementation, the constraint is NOT built into the DGA.
    % Part of this implementation is projecting the Jacobian onto the
    % constrained space.
    
    % Grab parameters.
    phi=q(1);
    theta=q(2);

    % Typical normal density assuming data is iid (w/in groups,
    % just independent between groups).
    dens = NormalDensity(cos(theta)*sin(phi),sin(theta)*sin(phi),...
        cos(phi),s1,s2,s3,dataX,dataY,dataZ);
    % The Jacobian calculation that projects onto the manifold
    jac = -log(JacobianParameterized(theta,phi,dataX,dataY,dataZ));
    % Note that I am omitting the Gram matrix here.  This density does not
    % assume the constraint at the level of the input.  I am not,
    % therefore, integrating over a function F that induces the manifold. I
    % am still, however, projecting the jacobian, since that derivative
    % should "know" that a constraint exists.
    nll = dens+jac;

    % The derivative of this density is a little tricky, so I'm doing it in
    % a helper function (below).
    dnll=get_dnll_density(q,s1,s2,s3,dataX,dataY,dataZ);
end



function dnll = get_dnll_density(q,s1,s2,s3,dataX,dataY,dataZ)
    % Method to get the derivative of the log likelihood.  I ended up
    % explicitedly calculating the derivative of the jacobian term (at the
    % end of this method).  All the commented stuff is the general matrix
    % algebra form of this derivative (my earlier, lazier, attempt).
    phi=q(1);
    theta=q(2);
    dnll=zeros(2,1);
    
    mu1=cos(theta)*sin(phi);
    mu2=sin(theta)*sin(phi);
    mu3=cos(phi);
    dTmu1=-sin(theta)*sin(phi);
    dTmu2=cos(theta)*sin(phi);
    dTmu3=0;
    dPmu1=cos(theta)*cos(phi);
    dPmu2=sin(theta)*cos(phi);
    dPmu3=-sin(phi);
    
    % First add in the values associated with the original density
    % function:
    for i=1:length(dataX)
        dnll(1)=dnll(1)-(dPmu1)*(1/s1^2)*(dataX(i)-mu1);
        dnll(1)=dnll(1)-(dPmu2)*(1/s2^2)*(dataY(i)-mu2);
        dnll(1)=dnll(1)-(dPmu3)*(1/s3^2)*(dataZ(i)-mu3);
    end
    for j=1:length(dataY)
        dnll(2)=dnll(2)-(dTmu1)*(1/s1^2)*(dataX(i)-mu1);
        dnll(2)=dnll(2)-(dTmu2)*(1/s2^2)*(dataY(i)-mu2);
        dnll(2)=dnll(2)-(dTmu3)*(1/s3^2)*(dataZ(i)-mu3);
    end
    
    
    % I am omitting the derivative of the Gram matrix.  Due to the method
    % of integration, here, I do not think that it is necessary.
    
    % Calculate the Jacobian matrix as one always does for the fiducial
    % approach
    % Full Jacobian matrix
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
    mat3(:,5)=dataZ-cos(phi)*ones(length(dataX),1);
    full_mat=[mat1;mat2;mat3];
    
    % Derivative wrt phi
    mat1 = zeros(length(dataX), 5);
    mat1(:,1)=(-cos(theta)*sin(phi)).*ones(length(dataX),1);
    mat1(:,2)=(-sin(theta)*cos(phi)).*ones(length(dataX),1);
    mat1(:,3)=-(cos(theta)*cos(phi)).*ones(length(dataX),1);
    mat2 = zeros(length(dataY), 5);
    mat2(:,1)=(-sin(theta)*sin(phi)).*ones(length(dataX),1);
    mat2(:,2)=(cos(theta)*cos(phi)).*ones(length(dataX),1);
    mat2(:,4)=-(sin(theta)*cos(phi)).*ones(length(dataX),1);
    mat3 = zeros(length(dataZ), 5);
    mat3(:,1)=(-cos(phi)).*ones(length(dataX),1);
    mat3(:,5)=sin(phi)*ones(length(dataX),1);
    dPfull_mat=[mat1;mat2;mat3];
    
    % Derivative wrt theta
    mat1 = zeros(length(dataX), 5);
    mat1(:,1)=(-sin(theta)*cos(phi)).*ones(length(dataX),1);
    mat1(:,2)=(-cos(theta)*sin(phi)).*ones(length(dataX),1);
    mat1(:,3)=(sin(theta)*sin(phi)).*ones(length(dataX),1);
    mat2 = zeros(length(dataY), 5);
    mat2(:,1)=(cos(theta)*cos(phi)).*ones(length(dataX),1);
    mat2(:,2)=(-sin(theta)*sin(phi)).*ones(length(dataX),1);
    mat2(:,4)=-(cos(theta)*sin(phi)).*ones(length(dataX),1);
    mat3 = zeros(length(dataZ), 5);
    dTfull_mat=[mat1;mat2;mat3];
    
    dnll(1)=dnll(1)-(1/2)*trace( ((full_mat')*full_mat)\(...
        (dPfull_mat')*full_mat +(full_mat')*dPfull_mat) );
    dnll(2)=dnll(2)-(1/2)*trace( ((full_mat')*full_mat)\(...
        (dTfull_mat')*full_mat +(full_mat')*dTfull_mat) );
end

