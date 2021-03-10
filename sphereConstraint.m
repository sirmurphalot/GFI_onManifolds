function [c,dc] = sphereConstraint(q, shift_vector)
    %Constraint, and derivative of the constraint, for the CHMC algorithm.
    x=q(1)-shift_vector(1);
    y=q(2)-shift_vector(2);
    z=q(3)-shift_vector(3);
    c=x^2+y^2+z^2-1;
    dc=[2*x;2*y;2*z]';
end