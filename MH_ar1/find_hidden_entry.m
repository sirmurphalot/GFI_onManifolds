function value = find_hidden_entry(sign_mat, skew_sym_mat)
%     sign_mat = 1 - sign_mat.*2;
    indices_to_grab = find(sign_mat==-1);
    if(isempty(indices_to_grab))
        value = -1;
    else
        if 1==prod(sign_mat)
            new_skew_sym_mat = skew_sym_mat(indices_to_grab,indices_to_grab);
            value = abs(det(new_skew_sym_mat));
        else
            value = -1;
        end
    end
    
end