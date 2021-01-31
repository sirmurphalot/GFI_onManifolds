function [nll,dnll] = nll_density(q,s1,s2,s3,dataX,dataY,dataZ)
    % Function to get the fiducial density for a (x,y) pair.
    % For this implementation, the constraint is NOT built into the DGA.
    % Part of this implementation is projecting the Jacobian onto the
    % constrained space.
    
    % Grab parameters.
    mu1=q(1);
    mu2=q(2);
    mu3=q(3);

    % Typical normal density assuming data is iid (w/in groups,
    % just independent between groups).
    dens = nl_NormalDensity(mu1,mu2,mu3,s1,s2,s3,dataX,dataY,dataZ);
    % The Jacobian calculation that projects onto the manifold
    jac = Jacobian(mu1,mu2,mu3,dataX,dataY,dataZ);
    % Note that I am omitting the Gram matrix here.  This density does not
    % assume the constraint at the level of the input.  I am not,
    % therefore, integrating over a function F that induces the manifold. I
    % am still, however, projecting the jacobian, since that derivative
    % should "know" that a constraint exists.
    nll = dens+jac;

    % The derivative of this density is a little tricky, so I'm doing it in
    % a helper function (below).
    dnll=get_dnll_density(q,s1,s2,s3,dataX,dataY,dataZ);
end



function dnll = get_dnll_density(q,s1,s2,s3,dataX,dataY,dataZ)
    % Method to get the derivative of the log likelihood.  I ended up
    % explicitedly calculating the derivative of the jacobian term (at the
    % end of this method).  All the commented stuff is the general matrix
    % algebra form of this derivative (my earlier, lazier, attempt).
    mu1=q(1);
    mu2=q(2);
    mu3=q(3);
    dnll=zeros(3,1);
    
    % First add in the values associated with the original density
    % function:
    for i=1:length(dataX)
        dnll(1)=dnll(1)-(1/s1^2)*(dataX(i)-mu1);
    end
    for j=1:length(dataY)
        dnll(2)=dnll(2)-(1/s2^2)*(dataY(j)-mu2);
    end
    for k=1:length(dataZ)
        dnll(3)=dnll(3)-(1/s3^2)*(dataZ(k)-mu3);
    end
    
    % I am omitting the derivative of the Gram matrix.  Due to the method
    % of integration, here, I do not think that it is necessary.
    
    % Finally, I add in the derivative of the Jacobian term:
    % Gather the matrices associated with the gradient of the implicit
    % function F (nabla_eta F).
%     nablaEtaF=[eye(2);[-mu1/sqrt(1-mu2^2-mu2^2) -mu2/sqrt(1-mu2^2-mu2^2)]];
%     dMu1NablaEtaF=[zeros(2,2);...
%         [-1/mu3 0]];
%     dMu2NablaEtaF=[zeros(2,2);...
%         [0 -1/mu3]];
%     dMu3NablaEtaF=[zeros(2,2);...
%         [mu1/mu3^2 mu2/mu3^2]];
    % I'm technically doing this derivative with the Chain Rule method
    % after inducing an implicit function F.  So, I have to divide off the
    % Gram matrix so that all that is left is my matrix U (which we proved
    % is equivalent to the Q matrix).  (...oof I may need to prove that 
    % the derivatives are the same as well.)
%     dMu3NegLogGramF=-1/mu3;
%     
%     % Gather the matrices associated with the gradient of the DGA G
%     % (nabla_theta G).
%     nablaThetaG=zeros(12,3);
%     nablaThetaG(1:4,1)=ones(4,1);
%     nablaThetaG(5:8,2)=ones(4,1);
%     nablaThetaG(9:12,3)=ones(4,1);
%     
    % Save space by defining the following dummy variable.
%     nablaGnablaF=nablaThetaG*nablaEtaF;
    
    % Add in the derivative information with respect to the jacobian term.
    % Note that there is no contribution to mu3.  This is ultimately a
    % consequence of nabla_theta G not being dependent on mu3 -- if this
    % were not the case, I believe we would have to include a term.
%     dnll(1)=dnll(1)+dMu1NegLogGramF;
%     dnll(1)=dnll(1)+trace((nablaGnablaF'*nablaGnablaF)\...
%         (nablaGnablaF'*nablaThetaG*dMu1NablaEtaF+...
%         dMu1NablaEtaF'*nablaThetaG'*nablaGnablaF));
% %     dnll(2)=dnll(2)+dMu2NegLogGramF;
%     dnll(2)=dnll(2)+trace((nablaGnablaF'*nablaGnablaF)\...
%         (nablaGnablaF'*nablaThetaG*dMu2NablaEtaF+...
%         dMu2NablaEtaF'*nablaThetaG'*nablaGnablaF));
%     dnll(3)=dnll(3)+dMu3NegLogGramF;
%     dnll(3)=dnll(3)+trace((nablaGnablaF'*nablaGnablaF)\...
%         (nablaGnablaF'*nablaThetaG*dMu3NablaEtaF+...
%         dMu3NablaEtaF'*nablaThetaG'*nablaGnablaF));

    % The above was the dirty way.  I've solved things out to get
    dnll(1)=dnll(1)-2*mu1/(mu1^2+mu2^2+mu3^2);
    dnll(2)=dnll(2)-2*mu2/(mu1^2+mu2^2+mu3^2);
    dnll(3)=dnll(3)+2*(mu1^2+mu2^2)/(mu3*(mu1^2+mu2^2+mu3^2));
%     dnll(3)=dnll(3)+2*mu3 - 1/mu3;
end

