function [c, ceq] = optim_constraint(a) %A_vector, Lambda
%     q = [A_vector';diag(Lambda)];
%     ceq = sum(abs(ar1_constraint(q)));
    c = -prod(1-2*a);
    
%     ceq =prod(1-2*a)-1;
    ceq = [];
    
    
end