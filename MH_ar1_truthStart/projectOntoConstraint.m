function [a, flag] = projectOntoConstraint(z, Tx, consFunc)
    %
    options = optimoptions('fsolve','Display','none','TolFun', ...
        8e-6, 'MaxFunctionEvaluations', 3e5, 'MaxIterations', 3e5);
    flag=1;
    solve_this = @(a) (consFunc(z + Tx*a));
    [~,d] = size(Tx);
    initial = zeros(d,1);
    a = fsolve(solve_this, initial, options);
    z_shift = z + Tx*a;
    
    tol=1e-2;
    if sum( abs(consFunc(z_shift)) < tol )==0
        flag = 0;
    end
end