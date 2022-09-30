function [y, prob, Tx, flag] = findProposal(curr_location, proposalScale, consFunc, dConsFunc)
% Author: Alexander Murph
    Qx = dConsFunc(curr_location);
    [Q,R] = qr(Qx);
    x = (R==0)';
    Tx = Q(:,x);
    [~,d] = size(Tx);
    center = zeros(d,1);
    scale = eye(d);
    
    v_temp = mvnrnd(center, proposalScale*scale, 1);
    v = Tx*v_temp';
    prob = mvnpdf(v_temp, center', proposalScale*scale);
    
    [a,flag] = projectOntoConstraint(curr_location+v, Qx, consFunc, []);
    y = curr_location + v + Qx*a;
end