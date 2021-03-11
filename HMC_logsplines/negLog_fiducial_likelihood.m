function nll = negLog_fiducial_likelihood(knots,coef_matrix,thetas,data)
    % Note that for now I'm assuming a vector of x data points. I might
    % generalize to multiple dimensions later, although I have yet to give
    % this much thought.
    [~,n] = size(data);
    
    nll = 0;
    for i=1:n
        nll_dens = nll_dens - ...
            log(logspline_density(data(1,i),knots,coef_matrix,thetas));
    end 
    nll = nll - Jacobian(knots,coef_matrix,thetas,data);  
end