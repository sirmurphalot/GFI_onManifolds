function [y, logdens, Nx, flag] = findProposal(curr_location, proposalScale, consFunc, dConsFunc)
    dc = dConsFunc(curr_location);
    [~,d] = size(dc);
%   P = eye(d) - ((dc')/(dc*(dc')))*dc;
%   [Q,R] = eig(P);
%   [Q1,R,~] = svds(P,2);
%   [Q,R,~] = svd(P);
    [Q,~] = qr(dc');

    % All singular values should be zeros or ones, since P is a projection 
    % matrix.  This corrects for numerical issues here.
    
    % Tx are the vectors orthogonal to the null space of G.  The vectors
    % that move OFF OF the tangent plane.
    %Nx = Q(:,~x);
    Nx = Q(:,1:(end-2));
    
    % Qx are the vectors parallel to the null space of G.  These vectors do
    % NOT move off of the tangent plane.
    Tx = Q(:,(end-1):end);
    
    [~,d] = size(Tx);
    center = zeros(d,1);
    scale = eye(d); 
 
    % I need to estimate the tangent vectors in the rho and sigma directions.
    % This is going to be dcurr_location/drho, dcurr_location/dsigma.
    delta = 1e-6;
    [curr_rho, curr_sigma] = AandLambdaToRhoandSigma(curr_location);
    dim=size(curr_location);
    
    [dA_dRho, dA_dSigma] = deriveParameterization(curr_rho, curr_sigma, curr_location);
    x=dim(1);
    n_timepoints = (1/2)*(-1+sqrt(1+8*x));

%     %temp_changeRho = RhoandSigmaToAandLambda(curr_rho+delta,curr_sigma,n_timepoints);
%     %A_changeRho = collapseParameters(temp_changeRho(1),temp_changeRho(2));
%     %temp_changeSigma = RhoandSigmaToAandLambda(curr_rho,curr_sigma+delta, n_timepoints);
%     %A_changeSigma = collapseParameters(temp_changeSigma(1),temp_changeSigma(2));
%   
% %     A_new = curr_location + (1-2.*rand([x,1]))*delta;
%     A_new = curr_location + (ones([x,1]))*delta;
%     [rho_new,sigma_new] = AandLambdaToRhoandSigma(A_new);
%  
  
%     dA_dSigma = dA_dSigma - ((dA_dSigma')*dA_dRho)./((dA_dRho')*dA_dRho).*dA_dRho;
% 
%     % Normalize these vectors.
    dA_dRho   = dA_dRho./(sqrt((dA_dRho')*dA_dRho));
    dA_dSigma = dA_dSigma./(sqrt((dA_dSigma')*dA_dSigma));
    dA_dRho   = dA_dRho - ((dA_dRho')*dA_dSigma)./((dA_dSigma')*dA_dSigma).*dA_dSigma;
%     dA_dSigma = dA_dSigma - ((dA_dSigma')*dA_dRho)./((dA_dRho')*dA_dRho).*dA_dRho;
    dA_dRho   = Tx*Tx'*(dA_dRho);
    dA_dSigma = Tx*Tx'*(dA_dSigma); 
     

    % Find the projections of each of these vectors onto the Tx space
    Tx_1 = Tx(:,1)./((Tx(:,1)')*Tx(:,1));
    Tx_2 = Tx(:,2)./((Tx(:,2)')*Tx(:,2));
    c11  = (dA_dRho)'*Tx_1;
    c12  = (dA_dRho)'*Tx_2;
    c21  = (dA_dSigma)'*Tx_1;
    c22  = (dA_dSigma)'*Tx_2;

    C_mat = [c11,c12;c21,c22];
%     scale(1,1) = 1;
    scale = real((C_mat')*scale*C_mat);% inv()
   % scale = eye(2);
    if(isnan(scale(1,1)))
        disp("hey!");
        flag = 0;
        logdens = 0;
        y = 0;
        return;
    end

    % Draw from a Cauchy distribution.
%     pd = makedist('tLocationScale','mu',0,'sigma',proposalScale,'nu',1);
%     v_temp = random(pd,1,d);
    v_temp = mvnrnd(center, proposalScale*scale, 1);
    
    % v should be a vector that moves ALONG the tangent plane.  Thus, we
    % make it a multiple of the vectors from Qx.
    v = Tx*v_temp';
%     logdens = sum(log(pdf(pd,v_temp)));
    dd = 2;
%     logdens = real(-(0.5*dd)*log(2*pi)-(0.5)*sum(log(eig(proposalScale*scale)))-(0.5)*((v_temp - center')/(proposalScale*scale))*(v_temp' - center));
     logdens = real(-(0.5*dd)*log(2*pi)-(0.5)*log(det(proposalScale*scale))-(0.5)*((v_temp - center')/(proposalScale*scale))*(v_temp' - center));
%      logdens = sum(log(mvnpdf(v_temp, center', proposalScale*scale)));
    
    % We now have a vector that is along the tangent plane.  We project it
    % downward using the vectors that go in the orthogonal directions.
%     [a,flag] = projectOntoConstraint(curr_location+v, Nx, consFunc, []);
%     
%     y = curr_location + v + Nx*a;

    [a,flag] = projectOntoConstraint(curr_location+v, dc', consFunc, []);
    
    y = curr_location + v + dc'*a;
    
    [tempA,tempL] = uncollapseParameters(y);
    temptempL = diag(tempL)
    if ~isequal(sort(diag(tempL)), diag(tempL)) && (flag)
         disp("check the proposal value");
         U = (eye(n_timepoints) - tempA)/(eye(n_timepoints) + tempA);
         sigma = U*(tempL.^2)*(U');
         [r,s] = AandLambdaToRhoandSigma(y)
         ar1_constraint(y)
    end
        
end
