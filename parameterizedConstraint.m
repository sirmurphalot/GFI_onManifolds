function [c,dc] = parameterizedConstraint(q)
    %Constraint, and derivative of the constraint, for the CHMC algorithm.
    phi=q(1);
    theta=q(2);
    if (phi<pi)&&(phi>0)&&(theta<=(2*pi))&&(theta>0)
        c=0;
    else
        c=Inf;
    end
    dc=[0;0]';
end