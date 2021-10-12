function [nll,dnll] = negLog_fiducial_likelihood(knots,thetas,data)
    % Note that for now I'm assuming a vector of x data points. I might
    % generalize to multiple dimensions later, although I have yet to give
    % this much thought.
    if sum(isnan(thetas))>0 || sum(isinf(thetas))>0 || sum(abs(thetas)>1e3)>0
       disp("getting bad theta values");
       nll = inf;
       dnll = zeros(size(thetas));
    else
        [n1,n2] = size(data);
        n = max(n1,n2);
        nll = 0;
        for i=1:n
            nll = nll - ...
                logspline_density(data(i),knots,thetas);
        end 
        [Jac, dJac] = Jacobian(knots,thetas,data);
        nll = nll - Jac;
        dnll = zeros(size(thetas));
        for k = 1:n
            for i=1:length(thetas)
                dnll(i) = -get_bspline_value(data(k),knots,i)-dJac(i);
            end
        end
    end
%     dnll = dnll';
end