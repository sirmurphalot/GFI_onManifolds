function [a, flag] = projectOntoConstraint(z, Tx, consFunc, curr_loc)
    %
    options = optimoptions('fsolve','Algorithm', 'levenberg-marquardt','Display',...
        'none', 'MaxFunctionEvaluations', 2e3, 'MaxIterations', 2e3);
    flag=1;
    solve_this = @(a) (consFunc(z + Tx*a));
    [~,d] = size(Tx);
    initial = zeros(d,1);
    a = fsolve(solve_this, initial, options);
    z_shift = z + Tx*a;
    tol=1e-2;
    
%     xx = length(z);
%     dd=(1/2)*(-1+sqrt(1+8*xx));
%     tol2 = zeros([1,xx]) + 0.1;
%     tol2((xx-dd+1):end) = tol2((xx-dd+1):end) + 0.4;
%     a_elements = z_shift(1:(xx-dd));
    
    if isempty(curr_loc)
        if sum( abs(consFunc(z_shift)) > tol ) > 0
            flag = 0;
            disp("projection failed!");
            disp(max(abs(consFunc(z_shift))));
            disp(sum( abs(consFunc(z_shift)) < tol ));
        end
        
    else
        if sum( (abs(z_shift-curr_loc)-tol) > 0 ) > 0 | sum( abs(consFunc(z_shift)) > tol ) > 0
            flag = 0;
            disp("reverse projection failed to get back to initial point!");
            disp(sum( abs(z_shift-curr_loc)  ));
            disp(sum( abs(consFunc(z_shift)) < tol ));
        end
    
    end
end
