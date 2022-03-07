function [c,ceq] = fminconstr(x)
    solve_this = @(a) (ar1_constraint(a));
    c = []; % No nonlinear inequality
    ceq = solve_this(x); % fsolve objective is fmincon nonlinear equality constraints
end