function [logprob, flag] = findReverseProposalProb(curr_loc, proposal, proposalScale, consFunc, dConsFunc)  
    % First calculate the projection matrix at the proposal.
    dc = dConsFunc(proposal);
    [~,d] = size(dc);
    P = eye(d) - ((dc')/(dc*(dc')))*dc;
    [Q,R,~] = svd(P);
    
    % All singular values should be zeros or ones, since P is a projection 
    % matrix.  This corrects for numerical issues here.
    x = (diag(R.^2)>=1e-4)';
    
    % Tx are the vectors orthogonal to the null space of G.  The vectors
    % that move OFF OF the tangent plane.
    Ty = Q(:,~x);
    % Qx are the vectors parallel to the null space of G.  These vectors do
    % NOT move off of the tangent plane.
    Qy = Q(:,x);
    
    % I want to move along the tangent space with my proposal, so I use Qy.
    [~,d] = size(Qy);
    center = zeros(d,1);
    scale = eye(d);
    
    x_minus_y = curr_loc - proposal;
    
    % Zero out the vector components that move OFF OF the tangent plane.
    % Note that the tangent vector that gets us above x on y's tangent plane is:
    % v' = x - y - w', where w' is orthogonal to y's tangent plane.
    v_temp = (Qy')*x_minus_y;
    pd = makedist('tLocationScale','mu',0,'sigma',proposalScale,'nu',1);
    logprob = sum(log(pdf(pd,v_temp')));
%     logprob = sum(log(mvnpdf(v_temp, center, proposalScale*scale)));

    % Now x_minus_y has been projected -- Qy Qy' = Py is a projection onto the
    % tangent plane at y.
    v_vector = Qy*v_temp;
    
    % Make sure that the return move is possible to find using our Newtong solver,
    % considering the curvature
    % of the manifold (it might not be, in which case we reject).
    [~, flag] = projectOntoConstraint(proposal + v_vector, Ty, consFunc);
end