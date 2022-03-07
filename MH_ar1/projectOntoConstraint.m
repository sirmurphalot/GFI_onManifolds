function [a, flag] = projectOntoConstraint(z, Tx, consFunc)
    %
    options = optimoptions('fsolve','Display','none','TolFun', ...
        7e-5, 'MaxFunctionEvaluations', 2e3, 'MaxIterations', 2e3);
    flag=1;
    solve_this = @(a) (consFunc(z + Tx*a));
    [~,d] = size(Tx);
    initial = zeros(d,1);
    a = fsolve(solve_this, initial, options);
    z_shift = z + Tx*a;
    
    tol=1e-3;
    if sum( abs(consFunc(z_shift)) < tol )==0
        flag = 0;
    end
end