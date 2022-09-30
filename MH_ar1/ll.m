function logDens = ll(rho, sigma, data_values)
    [n,d] = size(data_values);
    logDens = 0;
    starting_covariance = zeros(d, d);
    for i=1:d
        for j=1:d
            starting_covariance(i,j) = rho^(abs(i-j));
        end
    end
    starting_covariance = starting_covariance.*((sigma.^2)/(1-rho.^2));
    
    for i=1:n
        % Typical normal density assuming data is iid (w/in groups,
        % just independent between groups).
        logDens=logDens-(0.5*d)*log(2*pi)-(0.5)*sum(log(eig(starting_covariance)))-(0.5)*((data_values(i,:)')'/starting_covariance)*(data_values(i,:)'); 
    end

end