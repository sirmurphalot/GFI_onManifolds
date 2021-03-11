function N = NormalDensity(mu1,mu2,mu3,sd1,sd2,sd3,dataX,dataY,dataZ)
    % Calculate the density, compiled on the logscale the converted back,
    % of the given data.
    % This function is the same for both approaches: constrained DGA and
    % nonconstrained DGA.
    logDensity = 0;
    for i=1:length(dataX)
        logDensity = lognpdf(dataX(i), mu1, sd1) + logDensity;
    end
    for j=1:length(dataY)
        logDensity = lognpdf(dataY(j), mu2, sd2) + logDensity;
    end
    for k=1:length(dataZ)
        logDensity = lognpdf(dataZ(k), mu3, sd3) + logDensity;
    end
    N = -logDensity;
end