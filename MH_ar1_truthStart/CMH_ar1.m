% I'll start with simulating a data matrix.
rng(9)
n_timepoints = 10;
n_observations = 1;
data_matrix = zeros(n_observations, n_timepoints);
s_true = 0.9
rho = 0.8

for i=1:n_observations
    first_value = normrnd(0,s_true,1);
    data_matrix(i,1) = first_value;
    for j = 2:n_timepoints
        data_matrix(i,j) = rho*data_matrix(i,j-1) + normrnd(0,s_true,1);
    end
end

true_covariance = zeros(n_timepoints, n_timepoints);
for i=1:n_timepoints
    for j=1:n_timepoints
        true_covariance(i,j) = rho^(abs(i-j));
    end
end
true_covariance = true_covariance.*((s_true.^2)/(1-rho.^2))

solve_this = @(a) (ar1_constraint(a));
options = optimoptions('fsolve','Display','None','TolFun', ...
        8e-6, 'MaxFunctionEvaluations', 3e5, 'MaxIterations', 3e5);
[U, Lambda2, ~] = svd(true_covariance);
Lambda = diag(diag(Lambda2.^(0.5)));
original_sign = sign(det(U));
my_count = 0;
[d,~] = size(Lambda);
for i = 0:(2^d-1)
    temp_signs = de2bi(i,d);
    temp_signs = 2*temp_signs(1:d)-1;
    if(1==original_sign*prod(temp_signs))
        tempUhat = U*diag(temp_signs);
        tempA = (eye(d) - tempUhat)/(eye(d) + tempUhat);
        temp_vector = collapseParameters(tempA,Lambda);
%         initial = collapseParameters(actualA,Lambda);
        temp_vector = fsolve(solve_this, temp_vector, options);
        if( 1 == prod(temp_vector((d+1):end)>-1 & temp_vector((d+1):end)<1))
            disp("got here");
            [actualA,~] = uncollapseParameters(temp_vector);
            actualU = tempUhat;
            initial_guess = temp_vector
            break
        end
    end
end


%%

num_iters = 50000;
burn_in = 1;
proposal_scale = 5e-4

ll_function = @(x) ll_density(x, data_matrix);
constraintFunc = @(x) (ar1_constraint(x));
dConstraintFunc = @(x) (dar1_constraint(x));

warning('off');
[samples, accepts] = runConstrainedMH(num_iters, burn_in, ...
    proposal_scale, initial_guess, ll_function, ...
    constraintFunc, dConstraintFunc);
mean(accepts)
save("setseed_MH_ar1_Draws.mat",'samples');
save("MH_accepts_truevalue.mat", "accepts");
%% Use the samples to calculate values for rho and sigma2:

samples = load("setseed_MH_ar1_Draws.mat");
samples = samples.samples;
accepts = load("MH_accepts_truevalue.mat");
accepts = accepts.accepts;
mean(accepts);

[~,num_of_samples] = size(samples);
burnin = 2000;
rhos = zeros(1,num_of_samples-burnin);
s2s = zeros(1,num_of_samples-burnin); 

for i=1:(num_of_samples-burnin)
    [A, Lambda] = uncollapseParameters(samples(:,i+burnin));
    U = (eye(n_timepoints) - A)/(eye(n_timepoints) + A);
    sigma = U*(Lambda.^2)*(U');
    rhos(i) = sigma(1,2)/sigma(1,1);
    s2s(i) = sigma(1,1)*(1 - (sigma(1,2)/sigma(1,1)).^2);
end

histogram(rhos);
hold on
xline(rho, 'label', 'True Value', 'Color','r');
t = 'Samples of $\rho$ Parameter';
title(t,'interpreter','latex')
hold off

histogram(sqrt(s2s));
hold on
xline(s_true, 'label', 'True Value', 'Color', 'r');
t = 'Samples of $\sigma$ Parameter';
title(t,'interpreter','latex')
hold off