function [nll,dnll] = nll_density(q,s1,s2,s3,dataX,dataY,dataZ)
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
    jac = Jacobian(mu1,mu2,mu3,dataX,dataY,dataZ);
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

    % The above was the dirty way.  I've solved things out to get
    dnll(1)=dnll(1)-2*mu1/(mu1^2+mu2^2+mu3^2);
    dnll(2)=dnll(2)-2*mu2/(mu1^2+mu2^2+mu3^2);
    dnll(3)=dnll(3)+2*(mu1^2+mu2^2)/(mu3*(mu1^2+mu2^2+mu3^2));
end

