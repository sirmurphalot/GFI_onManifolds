function F = full_density_spherical_parameterized(phi,theta,...
    s1,s2,s3,dataX,dataY,dataZ)
    % Function to get the fiducial density for a (phi, theta) pair.
    % For this implementation, the constraint is built into the the DGA, so
    % we require no projection onto the constrained space.

    % To allow an array of phi, theta pairs to be used as an input,
    % we create a matrix of the correct size and iteratively collect
    % all density values.
    F = zeros(size(phi));
    for row=1:size(phi,1)
        for col=1:size(phi,2)
            % Typical normal density assuming data is iid (w/in groups,
            % just independent between groups).
            dens = nl_NormalDensity(sin(phi(row,col))*cos(theta(row,col)),...
                sin(phi(row,col))*sin(theta(row,col)),...
            cos(phi(row,col))+1,s1,s2,s3,dataX,dataY,dataZ);
            % The typical fiducial Jacobian without considering the
            % constraint
            jac = JacobianParameterized(theta(row,col),phi(row,col),...
                dataX,dataY,dataZ);
            F(row,col) = exp(-(dens+jac))*sin(phi(row,col));
        end
    end
end