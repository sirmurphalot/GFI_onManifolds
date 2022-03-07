function q = collapseParameters(A,Lambda)
    % Takes skew-symmetric matrix A and diagonal matrix Lambda and
    % collapses them into a single column vector of parameters q. 
    dim=size(Lambda);
    d=dim(1);
    q=zeros(d*(d+1)/2,1);
    
    % The first d(d-1)/2 parameters are of the A elements, read BY COLUMN
    % on the upper triangle.
    counter=1;
    for i=1:(d-1)
       for j=(i+1):d
           q(counter)=A(i,j);
           counter=counter+1;
       end
    end
    
    % The last d parameters in q are the diagonal elements of the Lambda
    % matrix.  
    q(counter:end)=diag(Lambda);
    
end