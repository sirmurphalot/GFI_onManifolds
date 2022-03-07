% env_id = "12"
rng(13)

n_timepoints = 50;
n_observations = 1;
data_matrix = zeros(n_observations, n_timepoints);
s_true = 0.9
rho = 0.7


first_value = normrnd(0,s_true/sqrt(1-rho.^2),1);
data_matrix(1) = first_value;
for j = 2:n_timepoints
    data_matrix(j) = rho*data_matrix(j-1) + normrnd(0,s_true,1);
end



% Conditional MLE estimates according to Hamilton.
sum_y_t = 0
sum_y_t_y_tm1 = 0
sum_y_tm1 = 0
sum_y_tm1_sq = 0

for i = 2:n_timepoints
    sum_y_t = sum_y_t + data_matrix(1,i);
    sum_y_t_y_tm1 = sum_y_t_y_tm1 + data_matrix(1,i-1)*data_matrix(1,i);
    sum_y_tm1 = sum_y_tm1 + data_matrix(1,i-1);
    sum_y_tm1_sq = sum_y_tm1_sq + data_matrix(1,i-1).^2;
end
data_vals_mat = [n_timepoints - 1, sum_y_tm1;sum_y_tm1, sum_y_tm1_sq];
data_vals_vec = [sum_y_t; sum_y_t_y_tm1];
estimates = data_vals_mat \ data_vals_vec;

rho_cond_mle = estimates(2,1);
sd2_cond_mle = 0
for i = 2:n_timepoints
    sd2_cond_mle = sd2_cond_mle + ((data_matrix(1,i) - estimates(1,1) - ...
        rho_cond_mle*data_matrix(1,i-1)).^2)/(n_timepoints-1);
end

starting_covariance = zeros(n_timepoints, n_timepoints);
for i=1:n_timepoints
    for j=1:n_timepoints
        starting_covariance(i,j) = rho_cond_mle^(abs(i-j));
    end
end
starting_covariance = starting_covariance.*((sd2_cond_mle)/(1-rho_cond_mle.^2));


options = optimoptions('fsolve','Display','None','TolFun', ...
        8e-6, 'MaxFunctionEvaluations', 3e5, 'MaxIterations', 3e5);
[U, Lambda2] = eig(starting_covariance);
origU = U;
Lambda = diag(diag(Lambda2.^(0.5)));
disp("About to calculate the A matrix.")
original_sign = sign(det(U));
my_count = 0;
count = 0;
[d,~] = size(Lambda);
found_starting_point = false;
for i = 0:(2^d-1)
    if mod( i , 5000 ) == 0
        disp(i);
    end
    temp_signs = de2bi(i,d);
    temp_signs = 2*temp_signs(1:d)-1;
    if(1==original_sign*prod(temp_signs))
        tempUhat = U*diag(temp_signs);
        
        E = eig(tempUhat);
        if ~(sum(abs(real(E)+1)<1e-4)==0)
            count = count + 1;
            continue;
        else
            disp("found one!")
            U = tempUhat;
            break
        end
    end
end

A = (eye(n_timepoints)-tempUhat)/(eye(n_timepoints)+tempUhat);

% Notice that 
starting_covariance(2,1)/starting_covariance(1,1)

% Yet
newU = (eye(n_timepoints)-A)/(eye(n_timepoints)+A);
newCov = newU * (Lambda.^2) * (newU');
newCov(2,1)/newCov(1,1)

% options = optimoptions('fsolve','Display','None','TolFun', ...
%         8e-6, 'MaxFunctionEvaluations', 3e5, 'MaxIterations', 3e5);
% [U, Lambda2] = eig(starting_covariance);
% Lambda = diag(diag(Lambda2.^(0.5)));
% disp("About to calculate the A matrix.")
% original_sign = sign(det(U));
% my_count = 0;
% count = 0;
% [d,~] = size(Lambda);
% found_starting_point = false;
% for i = 0:(2^d-1)
%     if mod( i , 5000 ) == 0
%         disp(i);
%     end
%     temp_signs = de2bi(i,d);
%     temp_signs = 2*temp_signs(1:d)-1;
%     if(1==original_sign*prod(temp_signs))
%         tempUhat = U*diag(temp_signs);
%         
%         E = eig(tempUhat);
%         if ~(sum(abs(real(E)+1)<1e-4)==0)
%             count = count + 1;
%             continue;
%         else
%             disp("found one!")
%             U = tempUhat;
%             break
%         end
%     end
% end
% 
% A = (eye(n_timepoints)+U)/(eye(n_timepoints)-U);
% 
% lb = [-1*ones(1,n_timepoints*(n_timepoints-1)/2),-1*ones(1,n_timepoints)*inf];
% ub = [ones(1,n_timepoints*(n_timepoints-1)/2),ones(1,n_timepoints)*inf];
% x0 = collapseParameters(A, Lambda);
% 
% disp("starting to find starting point numerically.")
% opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
% initial_guess = fmincon(@(x)0,x0,[],[],[],[],lb,ub,@fminconstr,opts);

% original_sign = sign(det(U));
% my_count = 0;
% count = 0;
% [d,~] = size(Lambda);
% found_starting_point = false;
% for i = 0:(2^d-1)
%     if mod( i , 5000 ) == 0
%         disp(i);
%     end
%     temp_signs = de2bi(i,d);
%     temp_signs = 2*temp_signs(1:d)-1;
%     if(1==original_sign*prod(temp_signs))
%         tempUhat = U*diag(temp_signs);
%         
%         E = eig(tempUhat);
%         if ~(sum(abs(real(E)+1)<1e-4)==0)
%             count = count + 1;
%             continue;
%         end
%         tempA = (eye(d) + tempUhat)/(eye(d) - tempUhat);
%         temp_vector = collapseParameters(tempA,Lambda);
%         
%         if( 1 == prod(temp_vector((d+1):end)>-1 & temp_vector((d+1):end)<1))
%             disp("got here");
%             [actualA,~] = uncollapseParameters(temp_vector);
%             actualU = tempUhat;
%             initial_guess = temp_vector
%             found_starting_point = true;
%             break
%         end
%     end
% end
% 
% if ~found_starting_point
%     error("couldn't find valid starting value")
% end
% 
% save("initial_guess.mat", "initial_guess");

% x = load("initial_guess.mat")
% initial_guess = x.initial_guess

num_iters = 5000;
burn_in = 1;
%%
proposal_scale = 5e-1

ll_function = @(x) ll_density(x, data_matrix);
constraintFunc = @(x) (ar1_constraint(x));
dConstraintFunc = @(x) (dar1_constraint(x));
disp("Entering the MCMC Algorithm.");
warning('off');
[samples, accepts] = runConstrainedMH(num_iters, burn_in, ...
    proposal_scale, initial_guess, ll_function, ...
    constraintFunc, dConstraintFunc);
mean(accepts)
save("setseed_MH_ar1_tighter_proj.mat",'samples');
save("MH_accepts_tighter_proj.mat", "accepts");
%% Use the samples to calculate values for rho and sigma2:

samples = load("setseed_MH_ar1_tighter_proj.mat");

accepts = load("MH_accepts_tighter_proj.mat");

% samples = load('setseed_MH_ar1_Draws.mat');
% accepts = load("MH_accepts_truevalue");
accepts = accepts.accepts;
samples = samples.samples;
mean(accepts);

[~,num_of_samples] = size(samples);
burnin = 1;
rhos = zeros(1,num_of_samples-burnin);
s2s = zeros(1,num_of_samples-burnin); 
s_diagonal_terms = zeros(1,num_of_samples-burnin); 

for i=1:(num_of_samples-burnin)
    [A, Lambda] = uncollapseParameters(samples(:,i+burnin));
%     [d,~] = size(Lambda);
    U = (eye(n_timepoints) - A)/(eye(n_timepoints) + A);
    sigma = U*(Lambda.^2)*(U');
%     rho_av = 0;
%     for j = 1:5
%         rho_av = rho_av + sigma(j,2)/sigma(j,1);
%     end
%     rhos(i) = rho_av/5;
    rhos(i) = sigma(1,2)/sigma(1,1);
    s2s(i) = sigma(1,1)*(1 - (sigma(1,2)/sigma(1,1)).^2);
    s_diagonal_terms(i) = sigma(1,1);
end

histogram(rhos);
hold on
xline(rho, 'label', 'True Value', 'Color','r');
xline(rho_cond_mle, 'label', 'MLE Estimate', 'Color','b');
t = 'Samples of $\rho$ Parameter';
title(t,'interpreter','latex')
hold off

histogram(sqrt(s2s));
hold on
xline(s_true, 'label', 'True Value', 'Color', 'r');
xline(sqrt(sd2_cond_mle), 'label', 'MLE Estimate', 'Color','b');
t = 'Samples of $\sigma$ Parameter';
title(t,'interpreter','latex')
hold off

histogram(s_diagonal_terms);
diagonal_true = s_true.^2/(1-rho.^2);
hold on
xline(diagonal_true, 'label', 'True Value', 'Color', 'r');
xline(sd2_cond_mle/(1-rho_cond_mle.^2), 'label', 'MLE Estimate', 'Color','b');
t = 'Samples of $\Sigma$ Parameter';
title(t,'interpreter','latex')
hold off