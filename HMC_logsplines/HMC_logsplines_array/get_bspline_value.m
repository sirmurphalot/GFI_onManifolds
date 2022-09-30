function val = get_bspline_value(y_val,knots,bspline_number)
    if (y_val >= knots(bspline_number)) && (y_val < knots(bspline_number+1))
        val = (y_val - knots(bspline_number))/(knots(bspline_number+1)-knots(bspline_number));
    elseif (y_val >= knots(bspline_number+1)) && (y_val < knots(bspline_number+2))
        val = (knots(bspline_number+2)-y_val)/(knots(bspline_number+2)-knots(bspline_number+1));
    else
        val = 0;
    end
end