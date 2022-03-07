function [y, logdens, Nx, flag] = findProposal(curr_location, proposalScale, consFunc, dConsFunc)
    dc = dConsFunc(curr_location);
    [~,d] = size(dc);
    P = eye(d) - ((dc')/(dc*(dc')))*dc;
    [Q,R,~] = svd(P);
    
    % All singular values should be zeros or ones, since P is a projection 
    % matrix.  This corrects for numerical issues here.
    x = (diag(R.^2)>=1e-4)';
    
    % Tx are the vectors orthogonal to the null space of G.  The vectors
    % that move OFF OF the tangent plane.
    Nx = Q(:,~x);
    % Qx are the vectors parallel to the null space of G.  These vectors do
    % NOT move off of the tangent plane.
    Tx = Q(:,x);
    
    [~,d] = size(Tx);
    center = zeros(d,1);
    scale = eye(d);
    
    % Draw from a Cauchy distribution.
    pd = makedist('tLocationScale','mu',0,'sigma',proposalScale,'nu',1);
    v_temp = random(pd,1,d);
%     v_temp = mvnrnd(center, proposalScale*scale, 1);
    
    % v should be a vector that moves ALONG the tangent plane.  Thus, we
    % make it a multiple of the vectors from Qx.
    v = Tx*v_temp';
    logdens = sum(log(pdf(pd,v_temp)));
%     logprob = sum(log(mvnpdf(v_temp, center', proposalScale*scale)));
    
    % We now have a vector that is along the tangent plane.  We project it
    % downward using the vectors that go in the orthogonal directions.
    
    
    
    % I'm trying to parameterized the projection process. 12/8/2021
    % I'll numerically solve for the values sd2, rho.  First, I'll set up
    % the parameters for the solver in matlab.
    [A, Lambda] = uncollapseParameters(curr_location);
    [n,~] = size(Lambda);
    solve_this = @(x)(parameterized_AR1_projection(x(1),x(2),Tx,curr_location,v,n));
    options = optimoptions('fsolve','Display','none','TolFun', ...
        1e-5, 'MaxFunctionEvaluations', 1e4, 'MaxIterations', 1e4);
    flag=1;

    % Next, I'll use the current location as the starting value for the
    % algorithm.
    U = (eye(n) - A)/(eye(n) + A);
    sigma = U*(Lambda.^2)*(U');
    rho = sigma(1,2)/sigma(1,1);
    s2s = sigma(1,1)*(1 - (sigma(1,2)/sigma(1,1)).^2);
    initial = [rho, s2s];

    % Now, I solve.
    [y,fval] = fsolve(solve_this, initial, options);

    if sum(fval)>1
%         flag = 0;
    end

    % With a solution, I can solve back for my covariance matrix, and in
    % turn get the proposed parameter vector.
    rho = y(1);
    sd2 = y(2);
    covariance_matrix = zeros(n, n);
    for i=1:n
        for j=1:n
            covariance_matrix(i,j) = rho^(abs(i-j));
        end
    end
    covariance_matrix = covariance_matrix.*((sd2)/(1-rho.^2));
    [U, Lambda2] = eig(covariance_matrix);
    Lambda = diag(diag(Lambda2.^(0.5)));

    %%% This is aligning the values from the sd2,rho.
    [Lambda_sorted, sorted_indices] = sort(diag(Lambda));
    U_sorted = U(:,sorted_indices);

    %%% This is aligning the values for the current location.
    [A_x, Lambda_x] = uncollapseParameters(curr_location);
    [~, sorted_indices_x] = sort(diag(Lambda_x));
    Ux = (eye(n) - A_x)/(eye(n) + A_x);
    Ux_sorted = Ux(:,sorted_indices_x);
%     Ax_sorted = (eye(n)+Ux_sorted)/(eye(n)-Ux_sorted);
%     curr_location_sorted = collapseParameters(Ax_sorted, diag(sorted_Lambda_x));

    %%% Now we determine a signature matrix for the sd2,rho value based on
    %%% its correlation with the eigenvecters for the current location.
%     [~,n_vectors] = size(Ux_sorted);
%     signature_matrix = zeros(n_vectors,1);
%     for k = 1:n_vectors
%         temp_sign = sign(corr(U_sorted(:,k), Ux_sorted(:,k)));
%         signature_matrix(k) = 1*temp_sign;
%     end
%     signature_matrix = diag(signature_matrix);
% 
%     U_final = U_sorted * signature_matrix;

    % I need to make sure that the eigenvalues are still all not negative
    % 1.
    U_final = U_sorted;
    E = eig(U_final);

    if ~(sum(abs(real(E)+1)<1e-5)==0)
        flag = 0;
    end

    % If they aren't, then I finish up :)
    I = eye(n);
    A = (I-U_final)/(I+U_final);
    y = collapseParameters(A,diag(Lambda_sorted));
end