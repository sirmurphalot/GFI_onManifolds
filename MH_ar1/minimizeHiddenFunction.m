function value = minimizeHiddenFunction(integer_value, skew_sym_matrix)
    [d,~]=size(skew_sym_matrix);
    temp_signs = de2bi(integer_value,d);
    temp_signs = 1-2*temp_signs(1:d);
    value = find_hidden_entry(temp_signs, skew_sym_mat);

end