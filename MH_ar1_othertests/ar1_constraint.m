function c = ar1_constraint(parameters)
    [A,Lambda] = uncollapseParameters(parameters);
    
    [d,~] = size(Lambda);
    U_mat = (eye(d)-A)/(eye(d)+A);
    cov_mat = U_mat*(Lambda.^2)*(U_mat');
    
    first_constraint_num = (d-1)*(d)/2;
    second_constraint_num = d-2;
    c = zeros(first_constraint_num + second_constraint_num,1);
    
%     disp("The rho val is currently:");
%     x = cov_mat(2,1)/cov_mat(1,1);
%     disp(x);
    
    count1 = 0;
    count2 = first_constraint_num;
    for j=1:(d-1)
        % Constraint that required the diagonal bands be equal.
        % These are recorded by COLUMN!
        % There are d(d-1)/2 of these.
        % They are put FIRST in constraint vector.
        for i=1:j
            count1 = count1 + 1;
            c(count1) = cov_mat(i,j) - cov_mat(i+1,j+1);
        end
        
        % Constraint that the top row proportion terms are the same.
        % These are recorded left to right, starting with the first
        % element.
        % There are d of these.
        % They are put SECOND in constraint vector.
        if(j~=(d-1))
            count2 = count2 + 1;
            c(count2) = cov_mat(1,j+1)/cov_mat(1,j) - cov_mat(1,j+2)/cov_mat(1,j+1);
        end
    end
end