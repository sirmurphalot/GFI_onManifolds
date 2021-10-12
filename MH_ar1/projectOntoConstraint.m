function [a, flag] = projectOntoConstraint(z, Qx, consFunc, dConsFunc)
%     options = optimoptions('fsolve','Display','iter-detailed','PlotFcn',@optimplotfirstorderopt,'TolFun', 1e-10);
    options = optimoptions('fsolve','Display','none','TolFun', 1e-10);
%     nmax=10000;
    tol=1e-6;
%     [~,d]=size(Qx);
%     a=zeros(d,1);
%     i=0;
    flag=1;
    solve_this = @(a) (consFunc(z + Qx*a));
    a = fsolve(solve_this, 0, options);
    z_shift = z + Qx*a;
    if norm(consFunc(z_shift)) > tol
        flag = 0;
    end
%     while norm(consFunc(z_shift)) > tol
% %         newton = @(deltaA) ((dConsFunc(z_shift)')*Qx*deltaA+...
% %             consFunc(z_shift)));
%         a = a - consFunc(z_shift)/((dConsFunc(z_shift)')*Qx);
% %         delta_a = fsolve(newton, delta_a_guess, options);
% %         a = a + delta_a;
% %         delta_a_guess = delta_a;
%         disp(norm(consFunc(z_shift)));
%         z_shift = z + Qx*a;
%         i = i+1;
%         if i > nmax
%             disp('got here');
%             flag = 0;
%             break
%         end
%     end
end