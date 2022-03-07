function [y, logdens, Nx, flag] = findProposal(curr_location, proposalScale, consFunc, dConsFunc)
    dc = dConsFunc(curr_location);
    [~,d] = size(dc);
    P = eye(d) - ((dc')/(dc*(dc')))*dc;
    [Q,R,~] = svd(P);
    
    % All singular values should be zeros or ones, since P is a projection 
    % matrix.  This corrects for numerical issues here.
    x = (diag(R.^2)>=1e-4)';
    
    % Tx are the vectors orthogonal to the null space of G.  The vectors
    % that move OFF OF the tangent plane.
    Nx = Q(:,~x);
    % Qx are the vectors parallel to the null space of G.  These vectors do
    % NOT move off of the tangent plane.
    Tx = Q(:,x);
    
    [~,d] = size(Tx);
    center = zeros(d,1);
    scale = eye(d);
    
    % Draw from a Cauchy distribution.
    pd = makedist('tLocationScale','mu',0,'sigma',proposalScale,'nu',1);
    v_temp = random(pd,1,d);
%     v_temp = mvnrnd(center, proposalScale*scale, 1);
    
    % v should be a vector that moves ALONG the tangent plane.  Thus, we
    % make it a multiple of the vectors from Qx.
    v = Tx*v_temp';
    logdens = sum(log(pdf(pd,v_temp)));
%     logprob = sum(log(mvnpdf(v_temp, center', proposalScale*scale)));
    
    % We now have a vector that is along the tangent plane.  We project it
    % downward using the vectors that go in the orthogonal directions.
    [a,flag] = projectOntoConstraint(curr_location+v, Nx, consFunc);
    
    y = curr_location + v + Nx*a;

end