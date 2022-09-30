function dc = dlogspline_constraint(knots, thetas)
    % Assuming c=0, we can cut back a bit on computations.
    [d1,d2] = size(thetas);
    d = max(d1,d2);
    % Assuming c=0, we can cut back a bit on computations.
    dc = zeros(d,1);
    for i=1:d
        dc(i) = bspline_expectation(knots,thetas,i);
    end
    dc = dc./sqrt(sum(dc.^2));
%     dc = dc';
end