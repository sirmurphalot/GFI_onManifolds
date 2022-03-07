% I'll start with simulating a data matrix.
env_id = getenv(char("SLURM_ARRAY_TASK_ID"))
rng(str2num(env_id))
upper_rho_output_file_name = strcat("UpperRhoCIs/UpperCI_arraynumber_", env_id, ".csv")
lower_rho_output_file_name = strcat("LowerRhoCIs/LowerCI_arraynumber_", env_id, ".csv")
upper_sd_output_file_name = strcat("UpperSdCIs/UpperCI_arraynumber_", env_id, ".csv")
lower_sd_output_file_name = strcat("LowerSdCIs/LowerCI_arraynumber_", env_id, ".csv")

n_timepoints = 20;
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
starting_covariance = starting_covariance.*((sd2_cond_mle.^2)/(1-rho_cond_mle.^2))

% solve_this = @(a) (ar1_constraint(a));
% options = optimoptions('fsolve','Display','None','TolFun', ...
%         8e-6, 'MaxFunctionEvaluations', 3e5, 'MaxIterations', 3e5);
[U, Lambda2] = eig(starting_covariance);
Lambda = diag(diag(Lambda2.^(0.5)));
original_sign = sign(det(U));
my_count = 0;
[d,~] = size(Lambda);
found_starting_point = false;
for i = 0:(2^d-1)
    temp_signs = de2bi(i,d);
    temp_signs = 2*temp_signs(1:d)-1;
    if(1==original_sign*prod(temp_signs))
        tempUhat = U*diag(temp_signs);
        % put an if statement here.
        tempA = (eye(d) - tempUhat)/(eye(d) + tempUhat);
        temp_vector = collapseParameters(tempA,Lambda);
        E = eig(tempUhat);
        if( 1 == prod(temp_vector((d+1):end)>-1 & temp_vector((d+1):end)<1) & sum(real(E)==-1)==0)
            disp("got here");
            [actualA,~] = uncollapseParameters(temp_vector);
            actualU = tempUhat;
            initial_guess = temp_vector
            found_starting_point = true;
            break
        end
    end
end

if ~found_starting_point
    error("couldn't find valid starting value")
end


%%

num_iters = 150;
burn_in = 1;
proposal_scale = 3e-6

ll_function = @(x) ll_density(x, data_matrix);
constraintFunc = @(x) (ar1_constraint(x));
dConstraintFunc = @(x) (dar1_constraint(x));

warning('off');
[samples, accepts] = runConstrainedMH(num_iters, burn_in, ...
    proposal_scale, initial_guess, ll_function, ...
    constraintFunc, dConstraintFunc);
% mean(accepts)
% save("setseed_MH_ar1_Draws.mat",'samples');
%% Use the samples to calculate values for rho and sigma2:

% samples = load("setseed_MH_ar1_Draws.mat").samples;
[~,num_of_samples] = size(samples);
burnin = 50;
rhos = zeros(1,num_of_samples-burnin);
s2s = zeros(1,num_of_samples-burnin); 

for i=1:(num_of_samples-burnin)
    [A, Lambda] = uncollapseParameters(samples(:,i+burnin));
    U = (eye(n_timepoints) - A)/(eye(n_timepoints) + A);
    sigma = U*(Lambda.^2)*(U');
    rhos(i) = sigma(1,2)/sigma(1,1);
    s2s(i) = sigma(1,1)*(1 - (sigma(1,2)/sigma(1,1)).^2);
end

proportion_rhos_above = mean(rhos>rho)
proportion_rhos_below = mean(rhos<rho)
proportion_sds_above = mean(sqrt(s2s)>s_true)
proportion_sds_below = mean(sqrt(s2s)<s_true)

writetable(table(proportion_rhos_below),upper_rho_output_file_name);

writetable(table(proportion_rhos_above),lower_rho_output_file_name);

writetable(table(proportion_sds_below),upper_sd_output_file_name);
writetable(table(proportion_sds_above),lower_sd_output_file_name);