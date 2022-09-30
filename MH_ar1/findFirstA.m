function difference = findFirstA(orthogonal_mat, A_vector, Lambda_vector)
    Lambda_vector = diag(Lambda_vector);
    q = [A_vector'; Lambda_vector];
    [A,Lambda] = uncollapseParameters(q);
    [~,d] = size(Lambda);
    cayley_transform = (eye(d)-A)/(eye(d)+A);
    mat_diff = orthogonal_mat - cayley_transform;
    mat_sq_diff = mat_diff.^2;
    difference = sqrt(sum(sum(mat_sq_diff)));
end