function G = full_density_spherical(phi,theta,s1,s2,s3,dataX,dataY,...
    dataZ, shift_vector)
    % Function to get the fiducial density for a (phi, theta) pair.
    % For this implementation, the constraint is NOT built into the DGA.
    % Part of this implementation is projecting the Jacobian onto the
    % constrained space.
    
    % To allow an array of (phi, theta) pairs to be used as an input,
    % we create a matrix of the correct size and iteratively collect
    % all density values.
    G = zeros(size(phi));
    for row=1:size(phi,1)
        for col=1:size(phi,2)
            % Typical normal density assuming data is iid (w/in groups,
            % just independent between groups).
            dens = nl_NormalDensity(sin(phi(row,col))*cos(theta(row,col))+shift_vector(1),...
                sin(phi(row,col))*sin(theta(row,col))+shift_vector(2),...
            cos(phi(row,col))+shift_vector(3),s1,s2,s3,dataX,dataY,dataZ);
            % The Jacobian calculation that projects onto the manifold
            jac = Jacobian(sin(phi(row,col))*cos(theta(row,col))+shift_vector(1),...
                sin(phi(row,col))*sin(theta(row,col))+shift_vector(2),...
                cos(phi(row,col))+shift_vector(3),...
                dataX,dataY,dataZ,shift_vector);
            % Note that I multiply by sin(phi).  This the the Gram matrix
            % multiplication that allows us to integrate on a submanifold.
            G(row,col) = exp( -( dens+jac ) )*sin(phi(row,col));
        end
    end
end