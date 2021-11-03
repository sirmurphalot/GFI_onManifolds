function [A,Lambda] = uncollapseParameters(q)
    % Takes skew-symmetric matrix A and diagonal matrix Lambda and
    % collapses them into a single column vector of parameters q. 
    
    % q should be of dimension (d(d+3))/2.  To find d, simply apply
    % the quadratic formula.
    dim=size(q);
    x=dim(1);
    d=(1/2)*(-1+sqrt(1+8*x));
    % Next, fill in A and Lambda in reverse of what is done in 
    % collapseParameters.m
    A=zeros(d,d);
    Lambda=diag(q(1:d));
    counter=1;
    for i=2:d
       for j=1:(i-1)
           A(i,j)=q(d+counter);
           A(j,i)=-q(d+counter);
           counter=counter+1;
       end
    end
end