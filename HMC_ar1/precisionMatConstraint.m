function [c,dc] = precisionMatConstraint(q)
    % I'll assume that this constraint is on a 4x4 matrix and the indicies
    % (2,4) are set to zero.
    zero_i=4;
    zero_j=2;

    % Start by gathering the parameters and instantiating the vectors:
    [A,Lambda]=uncollapseParameters(q);
    dim=size(A);
    IplusA = eye(dim(1))+A;
    IplusAinv=inv(IplusA);
    
    Derivatives=Derivator();
    Derivatives.dim=dim(1);
    Derivatives.Lambda=Lambda;
    Derivatives.IplusA=IplusA;
    Derivatives.IplusAinv=IplusAinv;
    dc=zeros(Derivatives.dim*(Derivatives.dim+1)/2,1);
    
    % We start by calculating the minor of the covariance matrix used to
    % get the general constraint c.
    CovarMat = ((IplusA')*IplusAinv)*(Lambda^2)*((IplusAinv')*IplusA);
    x=1:4;
    y=1:4;
    x(zero_i)=[];
    y(zero_j)=[];
    CovarMat_minor = CovarMat(x,y);

    % Fill in the derivatives with respect to Lambda.  There should be
    % dim many of them.
    for param=1:Derivatives.dim
        dc(param)=Derivatives.dLambda_g(zero_i,zero_j,param);
    end
    % Fill in the derivatives with respect to A.  There should be
    % dim(dim-1)/2 many of them.
    counter=1;
    for c=2:Derivatives.dim
        for d=1:(c-1)
            dc(Derivatives.dim+counter)=Derivatives.dA_g(zero_i,zero_j,c,d);
            counter=counter+1;
        end
    end
    
    % Finally, we return the relevant values.
    c = det(CovarMat_minor);
    dc=dc';
end