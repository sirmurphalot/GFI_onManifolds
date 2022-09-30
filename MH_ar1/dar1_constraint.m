function dc = dar1_constraint(parameters)
    % Assuming c=0, we can cut back a bit on computations.
    [A,Lambda] = uncollapseParameters(parameters);
    
    [d,~] = size(Lambda);
    U_mat = (eye(d)-A)/(eye(d)+A);
    cov_mat = (U_mat)*(Lambda.^2)*U_mat';
    IpA = eye(d) + A;
    IpAinv = inv(IpA);
    
    first_constraint_num = (d-1)*(d)/2;
    second_constraint_num = d-2;
    dc = zeros(first_constraint_num + second_constraint_num, length(parameters));
    
    
    %%% Calculate the A partials
    % The A partials are calculated BY ROW.
    countA = 0;
    for i=1:(d-1)
        for j=(i+1):d
            %%% Now, go through each constraint value in c
            count1 = 0;
            count2 = (d-1)*(d)/2;
            countA = countA + 1;
            
            %%% Calculate the matrix partial value
            Jij = zeros(d,d);
            Jij(i,j) = 1;
            B = IpAinv*(-Jij + Jij')*IpAinv*(Lambda.^2)*(IpAinv')*(IpA);
            cov_d_A_ij = 2*(B + B');
          
            %%% Fill in the dc vector with partials in the correct place.
            for l=1:(d-1)
                % First constraint, read by row.  There should be d(d-1)/2
                % of these.  
                for k=1:l
                    count1 = count1 + 1;
                    dc(count1, countA) = cov_d_A_ij(k,l) - cov_d_A_ij(k+1,l+1);
                end
                
                % Second constraint, read left to right.  There should be
                % d-1 of these.
                if(l~=(d-1))
                    count2 = count2 + 1;
                    %%% Calculate partial 
                    dc(count2, countA) = cov_mat(1,l+1)*cov_d_A_ij(1,l+1) + cov_d_A_ij(1,l+1)*cov_mat(1,l+1) - ...
                        cov_d_A_ij(1,l+2)*cov_mat(1,l) - cov_mat(1,l+2)*cov_d_A_ij(1,l);
                end
                
            end
            %%%%
        end
    end
    %%%
    
    
    %%% Calculate the Lambda partials
    countLambda = countA;
    for i=1:d
        %%% Now, go through each constraint value in c
        count1 = 0;
        count2 = (d-1)*(d)/2;
        countLambda = countLambda + 1;

        %%% Calculate the matrix partial value
        Jii = zeros(d,d);
        Jii(i,i) = 1;
        cov_d_Lambda_ii = 2*Lambda(i,i)*(IpA')*IpAinv*Jii*(IpAinv')*IpA;
        
        %%% Fill in the dc vector with partials in the correct place.
        for l=1:(d-1)
            
            % The first constraint, read by row.  There should be d(d-1)/2
            % of these.
            for k=1:l
                count1 = count1 + 1;
                dc(count1, countLambda) = cov_d_Lambda_ii(k,l) - cov_d_Lambda_ii(k+1,l+1);
            end
            
            % The second constraint, read left to right.  There should be
            % d-1 of these.
            if(l~=(d-1))
                count2 = count2 + 1;
                %%% Calculate partial 
                dc(count2, countLambda) = cov_mat(1,l+1)*cov_d_Lambda_ii(1,l+1) + cov_d_Lambda_ii(1,l+1)*cov_mat(1,l+1) - ...
                                            cov_d_Lambda_ii(1,l+2)*cov_mat(1,l) - cov_mat(1,l+2)*cov_d_Lambda_ii(1,l);
            end
            
        end
    end
    %%%%
    
end