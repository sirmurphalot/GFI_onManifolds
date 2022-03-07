function q = collapseParameters(A,Lambda)
    % Takes skew-symmetric matrix A and diagonal matrix Lambda and
    % collapses them into a single column vector of parameters q. 
    dim=size(Lambda);
    d=dim(1);
    q=zeros(d*(d+1)/2,1);
    
    % The first d parameters in q are the diagonal elements of the Lambda
    % matrix.  The next (d(d-1)/2) elements are from the lower triangle
    % of the matrix A.
    q(1:d,:)=diag(Lambda);
    counter=1;
    for i=2:d
       for j=1:(i-1)
           q(d+counter)=A(i,j);
           counter=counter+1;
       end
    end
end