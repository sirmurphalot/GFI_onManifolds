function [y, prob, Tx, flag] = findProposal(curr_location, proposalScale, consFunc, dConsFunc)
    dc = dConsFunc(curr_location);
    [~,d] = size(dc);
    P = eye(d) - ((dc')/(dc*(dc')))*dc;
    [Q,R,~] = svd(P);
    x = (diag(R.^2)>=1e-4)';
    
    % Tx are the vectors orthogonal to the null space of G.  The vectors
    % that move OFF OF the tangent plane.
    Tx = Q(:,~x);
    % Qx are the vectors parallel to the null space of G.  These vectors do
    % NOT move off of the tangent plane.
    Qx = Q(:,x);
    [~,d] = size(Qx);
    center = zeros(d,1);
    scale = eye(d);
    
    v_temp = mvnrnd(center, proposalScale*scale, 1);
    
    % v should be a vector that moves ALONG the tangent plane.  Thus, we
    % make it a multiple of the vectors from Qx.
    v = Qx*v_temp';
    prob = mvnpdf(v_temp, center', proposalScale*scale);
    
    % We now have a vector that is along the tangent plane.  We project it
    % downward using the vectors that go in the orthogonal directions.
    [a,flag] = projectOntoConstraint(curr_location+v, Tx, consFunc);
    
    y = curr_location + v + Tx*a;
end