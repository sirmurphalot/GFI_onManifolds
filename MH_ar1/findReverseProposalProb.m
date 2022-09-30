function [logprob, flag] = findReverseProposalProb(curr_loc, proposal, proposalScale, consFunc, dConsFunc)  
    % First calculate the projection matrix at the proposal.
    dc = dConsFunc(proposal);
    [other_d,d] = size(dc);
    P = eye(d) - ((dc')/(dc*(dc')))*dc;
%     [Q1,R,~] = svds(P,2);
%     [Q,R,~] = svd(P);
    %[Q,R] = eig(P,'nobalance');   
    [Q,~] = qr(dc');

    % Tx are the vectors orthogonal to the null space of G.  The vectors
    % that move OFF OF the tangent plane.
    Ty = Q(:,1:(end-2));
    
    % Qx are the vectors parallel to the null space of G.  These vectors do
    % NOT move off of the tangent plane.
    Qy = Q(:,(end-1):end);
    
    % I want to move along the tangent space with my proposal, so I use Qy.
    [~,d] = size(Qy);
    center = zeros(d,1);
    scale = eye(d);
    % This is where I force the process to move further in the rho direction:
%     scale(1,1) = 2;

    % I need to estimate the tangent vectors in the rho and sigma directions.
    % This is going to be dcurr_location/drho, dcurr_location/dsigma.
    delta = 1e-6;
    [curr_rho, curr_sigma] = AandLambdaToRhoandSigma(proposal);
    [dA_dRho, dA_dSigma] = deriveParameterization(curr_rho, curr_sigma, proposal);
    dim=size(proposal);
    x=dim(1);
    n_timepoints = (1/2)*(-1+sqrt(1+8*x));
    
%     % Normalize these vectors.
    dA_dRho   = dA_dRho./(sqrt((dA_dRho')*dA_dRho));
    dA_dSigma = dA_dSigma./(sqrt((dA_dSigma')*dA_dSigma));
    dA_dRho   = dA_dRho - ((dA_dRho')*dA_dSigma)./((dA_dSigma')*dA_dSigma).*dA_dSigma;
%     dA_dSigma = dA_dSigma - ((dA_dSigma')*dA_dRho)./((dA_dRho')*dA_dRho).*dA_dRho;
    dA_dRho   = Qy*Qy'*(dA_dRho);
    dA_dSigma = Qy*Qy'*(dA_dSigma); 
    

    % Find the projections of each of these vectors onto the Tx space
    Tx_1 = Qy(:,1)./((Qy(:,1)')*Qy(:,1));
    Tx_2 = Qy(:,2)./((Qy(:,2)')*Qy(:,2));
    c11  = (dA_dRho)'*Tx_1;
    c12  = (dA_dRho)'*Tx_2;
    c21  = (dA_dSigma)'*Tx_1;
    c22  = (dA_dSigma)'*Tx_2;

    T_proposal_y = (Qy/(Qy'*Qy))*Qy';
    
    C_mat = [c11,c12;c21,c22];
%     scale(1,1) = 1;
    scale = real((C_mat')*scale*C_mat); %inv()
    if(isnan(scale(1,1)))
        disp("hey")
        logprob = 0;
        flag = 0;
        y = 0;
        return;
    end
    %scale = eye(2);

    % Now, back to the original algorithm.   
    x_minus_y = curr_loc - proposal;
    
    % Zero out the vector components that move OFF OF the tangent plane.
    % Note that the tangent vector that gets us above x on y's tangent plane is:
    % v' = x - y - w', where w' is orthogonal to y's tangent plane.
    v_temp = (Qy')*x_minus_y;
%     pd = makedist('tLocationScale','mu',0,'sigma',proposalScale,'nu',1);
%     logprob = sum(log(pdf(pd,v_temp')));

    dd = 2;
    logprob = real(-(0.5*dd)*log(2*pi)-(0.5)*log(det(proposalScale*scale))-(0.5)*((v_temp' - center')/(proposalScale*scale))*(v_temp - center));
%      logprob = sum(log(mvnpdf(v_temp, center, proposalScale*scale)));

%     if logprob < -50
%        disp("wait!!"); 
%     end

    % Now x_minus_y has been projected -- Qy Qy' = Py is a projection onto the
    % tangent plane at y.
    v_vector = Qy*v_temp;
    
    % Make sure that the return move is possible to find using our Newtong solver,
    % considering the curvature
    % of the manifold (it might not be, in which case we reject).
%     [~, flag] = projectOntoConstraint(proposal + v_vector, Ty, consFunc, curr_loc);
    [~, flag] = projectOntoConstraint(proposal + v_vector, dc', consFunc, curr_loc);
end
