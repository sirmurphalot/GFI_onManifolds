function [prob, flag] = findReverseProposalProb(curr_loc, proposal, proposalScale, Tx, consFunc, dConsFunc)
% Author: Alexander Murph
    tol=1e-6;
    options = optimoptions('fsolve','Display','none','TolFun', 1e-10);
    
    Qy = dConsFunc(proposal);
    [Q,R] = qr(Qy);
    cols = (R==0)';
    Ty = Q(:,cols);
    [~,d] = size(Ty);
    center = zeros(d,1);
    scale = eye(d);
    
    x_minus_y = curr_loc - proposal;
    v_temp = (Ty')*x_minus_y;
    prob = mvnpdf(v_temp', center', proposalScale*scale);
    
    v_vector = Ty*v_temp;
    solve_this = @(a) (consFunc(proposal + v_vector + Qy*a));
    try
        a = fsolve(solve_this, 0, options);
    catch
        disp("hey");
    end
    
    z_shift = proposal + v_vector + Qy*a;
    
    flag=1;
    if norm(z_shift - curr_loc) > tol
        flag = 0;
    end
    
    
    if prob == 0
%         disp("hey");
    end
    
    % Also check inequality constraint:
%     h_func = det(Ty'*Tx);
%     if h_func < 0
%         flag = 0;
%     end
    
end