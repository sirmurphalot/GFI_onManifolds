function [nll,dnll] = nll_density(q,s1,s2,s3,dataX,dataY,dataZ,shift_vector)
    % Function to get the fiducial density for a (x,y) pair.
    % For this implementation, the constraint is NOT built into the DGA.
    % Part of this implementation is projecting the Jacobian onto the
    % constrained space.
    
    % Grab parameters.
    mu1=q(1);
    mu2=q(2);
    mu3=q(3);

    % Typical normal density assuming data is iid (w/in groups,
    % just independent between groups).
    dens = NormalDensity(mu1,mu2,mu3,s1,s2,s3,dataX,dataY,dataZ);
    % The Jacobian calculation that projects onto the manifold
    jac = Jacobian(mu1,mu2,mu3,dataX,dataY,dataZ,shift_vector);
    % Note that I am omitting the Gram matrix here.  This density does not
    % assume the constraint at the level of the input.  I am not,
    % therefore, integrating over a function F that induces the manifold. I
    % am still, however, projecting the jacobian, since that derivative
    % should "know" that a constraint exists.
    nll = dens+jac; %+log(mu3-shift_vector(3))

    % The derivative of this density is a little tricky, so I'm doing it in
    % a helper function (below).
    dnll=get_dnll_density(q,s1,s2,s3,dataX,dataY,dataZ,shift_vector);
end



function dnll = get_dnll_density(q,s1,s2,s3,dataX,dataY,dataZ,shift_vector)
    mu1=q(1);
    mu2=q(2);
    mu3=q(3);
    dnll=zeros(3,1);
    
    % First add in the values associated with the original density
    % function:
    for i=1:length(dataX)
        dnll(1)=dnll(1)-(1/s1^2)*(dataX(i)-mu1);
    end
    for j=1:length(dataY)
        dnll(2)=dnll(2)-(1/s2^2)*(dataY(j)-mu2);
    end
    for k=1:length(dataZ)
        dnll(3)=dnll(3)-(1/s3^2)*(dataZ(k)-mu3);
    end
    
    % I am omitting the derivative of the Gram matrix.  Due to the method
    % of integration, here, I do not think that it is necessary.
    mu1 = mu1 - shift_vector(1);
    mu2 = mu2 - shift_vector(2);
    mu3 = mu3 - shift_vector(3);
    P_mat = [mu2^2 + mu3^2, -mu1*mu2, -mu1*mu3;...
        -mu1*mu2, mu1^2+mu3^2, -mu2*mu3;...
        -mu1*mu3, -mu2*mu3, mu1^2+mu2^2];

    dm1P = [0, -mu2, -mu3;-mu2, 2*mu1, 0;-mu3, 0, 2*mu1];
    dm2P = [2*mu2, -mu1, 0;-mu1, 0, -mu3;0, -mu3, 2*mu2];
    dm3P = [2*mu3, 0, -mu1;0, 2*mu3, -mu2;-mu1, -mu2, 0];
    
    mat1 = zeros(length(dataX), 3);
    mat1(:,1)=ones(length(dataX),1);
    mat2 = zeros(length(dataY), 3);
    mat2(:,2)=ones(length(dataY),1);
    mat3 = zeros(length(dataZ), 3);
    mat3(:,3)=ones(length(dataZ),1);
    Jac_mat=[mat1;mat2;mat3];
    penrose_inverse = pinv(Jac_mat*P_mat*(Jac_mat'));
   
    dnll(1)=dnll(1)-trace(penrose_inverse*Jac_mat*dm1P*(Jac_mat'));
    dnll(2)=dnll(2)-trace(penrose_inverse*Jac_mat*dm2P*(Jac_mat'));
    dnll(3)=dnll(3)-trace(penrose_inverse*Jac_mat*dm3P*(Jac_mat'));
    % The above was the dirty way.  I've solved things out to get
%     dnll(1)=dnll(1)-2*mu1/(mu1^2+mu2^2+mu3^2);
%     dnll(2)=dnll(2)-2*mu2/(mu1^2+mu2^2+mu3^2);
%     dnll(3)=dnll(3)+2*(mu1^2+mu2^2)/(mu3*(mu1^2+mu2^2+mu3^2));
    
end

