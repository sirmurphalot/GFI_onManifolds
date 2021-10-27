function [a, flag] = projectOntoConstraint(z, Tx, consFunc)
    options = optimoptions('fsolve','Display','none','TolFun', ...
        1e-9, 'MaxFunctionEvaluations', 1e5, 'MaxIterations', 1e5);
    flag=1;
    solve_this = @(a) (consFunc(z + Tx*a));
    [~,d] = size(Tx);
    initial = zeros(d,1);
    a = fsolve(solve_this, initial, options);
    z_shift = z + Tx*a;
    
    tol=1e-1;
    if norm(consFunc(z_shift)) > tol
        flag = 0;
    end
end