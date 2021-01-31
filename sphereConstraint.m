function [c,dc] = sphereConstraint(q)
    %Constraint, and derivative of the constraint, for the CHMC algorithm.
    x=q(1);
    y=q(2);
    z=q(3)-1;
    c=x^2+y^2+z^2-1;
    dc=[2*x;2*y;2*z]';
end