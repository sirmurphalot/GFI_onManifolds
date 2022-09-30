function [dLoc_drho, dLoc_dsigma] = deriveParameterization(curr_rho, curr_sigma, curr_location)
    [n,~] = size(curr_location);
    dLoc_drho = zeros(n,1);
    dLoc_dsigma = zeros(n,1);
    d=(1/2)*(-1+sqrt(1+8*n));
    [A,Lambda] = uncollapseParameters(curr_location);
    IpAinv = inv(eye(d) + A);
    IpA = eye(d) + A;
    
    count = 1;
    for i = 1:d
        for j = (i+1):d
            % Get the A partials
            Jij = zeros(d,d);
            Jij(i,j) = 1;
            B = IpAinv*(-Jij + Jij')*IpAinv*(Lambda.^2)*(IpAinv')*(IpA);
            cov_d_A_ij = 2*(B + B');
            
            total_partial_wrt_rho = 0;
            total_partial_wrt_sigma = 0;
            for ii = 1:d
                for jj = i:d
                    dSigma_dRho = curr_sigma^2 * (curr_rho^(jj-ii) + (1-(jj-ii))*curr_rho^(jj-ii))/((1-curr_rho)^2) ;
                    total_partial_wrt_rho = total_partial_wrt_rho + ...
                        dSigma_dRho * cov_d_A_ij(ii,jj)^(-1);
                    dSigma_dsigma = 2*curr_sigma * curr_rho^(jj-ii) / (1-curr_rho);
                    total_partial_wrt_sigma = total_partial_wrt_sigma + ...
                        dSigma_dsigma * cov_d_A_ij(ii,jj)^(-1);
                end
            end
            dLoc_drho(count) = total_partial_wrt_rho;
            dLoc_dsigma(count) = total_partial_wrt_sigma;
            count = count + 1;
        end
    end
        
    % Get the Lambda partials
    for i = 1:d  
        Jii = zeros(d,d);
        Jii(i,i) = 1;
        cov_d_Lambda_i = 2*Lambda(i,i)*(IpA')*IpAinv*Jii*(IpAinv')*IpA;
        
        total_partial_wrt_rho = 0;
        total_partial_wrt_sigma = 0;
        for ii = 1:d
            for jj = i:d
                dSigma_dRho = curr_sigma^2 * (curr_rho^(jj-ii) + (1-(jj-ii))*curr_rho^(jj-ii))/((1-curr_rho)^2) ;
                total_partial_wrt_rho = total_partial_wrt_rho + ...
                    dSigma_dRho * cov_d_Lambda_i(ii,jj)^(-1);
                dSigma_dsigma = 2*curr_sigma * curr_rho^(jj-ii) / (1-curr_rho);
                total_partial_wrt_sigma = total_partial_wrt_sigma + ...
                    dSigma_dsigma * cov_d_Lambda_i(ii,jj)^(-1);
            end
        end
        dLoc_drho(count) = total_partial_wrt_rho;
        dLoc_dsigma(count) = total_partial_wrt_sigma;
        count = count + 1;
    end
end