function [nll,dnll] = negLog_fiducial_likelihood(knots,coef_matrix,thetas,data)
    % Note that for now I'm assuming a vector of x data points. I might
    % generalize to multiple dimensions later, although I have yet to give
    % this much thought.
    [~,n] = size(data);
    
    nll = 0;
    for i=1:n
        nll = nll - ...
            log(logspline_density(data(1,i),knots,coef_matrix,thetas));
    end 
    [Jac, dJac] = Jacobian(knots,coef_matrix,thetas,data);
    nll = nll - Jac;
    dnll = zeros(size(thetas));
    for i=1:length(thetas)
        dnll(i) = -sum(get_bspline_value(data,knots,coef_matrix,i))-dJac(i);
    end
%     dnll = dnll';
end