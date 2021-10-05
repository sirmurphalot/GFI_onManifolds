function [nll,dnll] = negLog_fiducial_likelihood(knots,thetas,data)
    % Note that for now I'm assuming a vector of x data points. I might
    % generalize to multiple dimensions later, although I have yet to give
    % this much thought.
    [~,n] = size(data);
    
    nll = 0;
    for i=1:n
        nll = nll - ...
            logspline_density(data(i),knots,thetas);
    end 
    [Jac, dJac] = Jacobian(knots,thetas,data);
    nll = nll - Jac;
    dnll = zeros(size(thetas));
    for i=1:length(thetas)
        dnll(i) = -get_bspline_value(data,knots,i)-dJac(i);
    end
%     dnll = dnll';
end