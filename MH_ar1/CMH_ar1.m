% I'll start with simulating a data matrix.
n_timepoints = 5;
n_observations = 50;
data_matrix = zeros(n_observations, n_timepoints);
sigma2 = 0.5;
rho = 0.1;

for i=1:n_observations
    first_value = normrnd(0,sqrt(sigma2),1);
    data_matrix(i,1) = first_value;
    for j = 2:n_timepoints
        data_matrix(i,j) = rho*data_matrix(i,j-1) + normrnd(0,sqrt(sigma2),1);
    end
end


[U, Lambda, ~] = svd(cov(data_matrix));

original_sign = sign(det(U));
my_count = 0;
[d,~] = size(Lambda);
for i = 0:(2^d-1)
    temp_signs = de2bi(i,8);
    temp_signs = 2*temp_signs(1:d)-1;
    if(1==original_sign*prod(temp_signs))
        tempUhat = U*diag(temp_signs);
        tempA = (eye(d) - tempUhat)/(eye(d) + tempUhat);
        fake_Lambda = zeros(d);
        temp_vector = collapseParameters(tempA,fake_Lambda);
        if(1 == prod(temp_vector>-1 & temp_vector<1))
            disp("got here");
            [actualA,~] = uncollapseParameters(temp_vector);
            actualU = tempUhat;
        end
    end
end


solve_this = @(a) (ar1_constraint(a));
options = optimoptions('fsolve','Display','none','TolFun', 1e-10);
initial = collapseParameters(actualA,Lambda);
% initial = [-5, log(0.1), log(0.4), log(1.6), log(1.5), log(2.3), log(1.5), log(2.4), log(2.1),...
%     log(1.2), log(1.25), log(1.3), log(0.6), log(0.7),...
%     log(0.5), log(0.4), log(0.3), ...
%     -2.5];
% initial = [-2, log(1), log(2), log(1), log(0.9), log(0.7), log(0.4), -2];
initial_guess = fsolve(solve_this, initial, options)

num_iters = 10000;
burn_in = 1;
proposal_scale = 3e-6

ll_function = @(x) ll_density(x, data_matrix);
constraintFunc = @(x) (ar1_constraint(x));
dConstraintFunc = @(x) (dar1_constraint(x));

warning('off');
[samples, accepts] = runConstrainedMH(num_iters, burn_in, ...
    proposal_scale, initial_guess, ll_function, ...
    constraintFunc, dConstraintFunc);
mean(accepts)
save("MH_ar1_Draws.mat",'samples');
%% Use the samples to calculate values for rho and sigma2:

samples = load("MH_ar1_Draws.mat").samples;
[~,num_of_samples] = size(samples);
rhos = zeros(1,num_of_samples);
s2s = zeros(1,num_of_samples); 

for i=1:num_of_samples
    [A, Lambda] = uncollapseParameters(samples(:,i));
    U = (eye(n_timepoints) - A)/(eye(n_timepoints) + A);
    sigma = U*(Lambda.^2)*(U');
    rhos(i) = sigma(1,2)/sigma(1,1);
    s2s(i) = sigma(1,1)/(1 - (sigma(1,2)/sigma(1,1)).^2);
end

histogram(rhos);
hold on
xline(rho);
hold off

histogram(s2s);
hold on
xline(sigma2);
hold off