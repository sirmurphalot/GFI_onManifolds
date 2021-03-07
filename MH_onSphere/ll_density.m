function ll = ll_density(q,s1,s2,s3,dataX,dataY,dataZ,z_shift)
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
    jac = Jacobian(mu1,mu2,mu3,dataX,dataY,dataZ,z_shift);
    % Note that I am omitting the Gram matrix here.  This density does not
    % assume the constraint at the level of the input.  I am not,
    % therefore, integrating over a function F that induces the manifold. I
    % am still, however, projecting the jacobian, since that derivative
    % should "know" that a constraint exists.
    ll = dens+jac - log(mu3);
end



