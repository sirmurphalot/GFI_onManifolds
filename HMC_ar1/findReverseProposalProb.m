function [prob, flag] = findReverseProposalProb(curr_loc, proposal, proposalScale, consFunc, dConsFunc)    
    dc = dConsFunc(proposal);
    [~,d] = size(dc);
    P = eye(d) - ((dc')/(dc*(dc')))*dc;
    [Q,R,~] = svd(P);
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
    v_temp = (Qy')*x_minus_y;
    pd = makedist('tLocationScale','mu',0,'sigma',proposalScale,'nu',1);
    prob = log(prod(pdf(pd,v_temp)));
%     prob = mvnpdf(v_temp, center, proposalScale*scale);
    v_vector = Qy*v_temp;
    
    if prob == 0
        disp("hey");
    end
    
    % Make sure that the return move is possible, considering the curvature
    % of the manifold (it might not be, in which case we reject).
    [a, flag] = projectOntoConstraint(proposal + v_vector, Ty, consFunc);
    z_shift = proposal + v_vector + Ty*a;
    
    % I need to know that the return projection got back to the original
    % value.
    tol=1e-1;
%     if norm(z_shift-curr_loc) > tol
%         disp("return probability failed!");
%         flag = 0;
%     end
    
end