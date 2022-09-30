function [rho,sigma] = AandLambdaToRhoandSigma(q)
    [A,Lambda] = uncollapseParameters(q);
    [n_timepoints,~] = size(A);
    U = (eye(n_timepoints) - A)/(eye(n_timepoints) + A);
    sigma = U*(Lambda.^2)*(U');
    rho = sigma(1,2)/sigma(1,1);
    s2s  = sigma(1,1)*(1 - (sigma(1,2)/sigma(1,1)).^2);
    sigma = sqrt(s2s);
end
