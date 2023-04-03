function [a, flag] = projectOntoConstraint(z, Tx, consFunc, curr_loc)
    %
    options = optimoptions('fsolve','Display','none','TolFun', ...
        1e-7, 'MaxFunctionEvaluations', 1e4, 'MaxIterations', 1e4);
    flag=1;
    solve_this = @(a) (consFunc(z + Tx*a));
    [~,d] = size(Tx);
    initial = zeros(d,1);
    a = fsolve(solve_this, initial, options);
    z_shift = z + Tx*a;
    tol=5e-1;
    
    xx = length(z);
    dd=(1/2)*(-1+sqrt(1+8*xx));
    
    if isempty(curr_loc)
        if sum( abs(consFunc(z_shift)) > tol ) > 0
            flag = 0;
            disp("projection failed!");
            disp(max(abs(consFunc(z_shift))));
            disp(sum( abs(consFunc(z_shift)) < tol ));
        end
        
    else
        if sum( (sum(abs(z_shift-curr_loc))-tol) > 0 ) > 0 | sum( abs(consFunc(z_shift)) > tol ) > 0
            flag = 0;
            disp("reverse projection failed to get back to initial point!");
            disp(sum( abs(z_shift-curr_loc)  ));
            disp(sum( abs(consFunc(z_shift)) < tol ));
        end
    
    end
end
