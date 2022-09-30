env_id = getenv(char("SLURM_ARRAY_TASK_ID"))
% env_id = "1"
rng(str2num(env_id))

n_timepoints = 30;
n_observations = 1;
data_matrix = zeros(n_observations, n_timepoints);
s_true = 1
rho = 0.7
sigma_true = s_true.^2/1-rho.^2


my_file_name_data = strcat("data_logs/sample_data_", env_id, ".mat");


samples = load(my_file_name_data);
samples = samples.samples;
%% Use the samples to calculate values for rho and sigma2:


[~,num_of_samples] = size(samples);
burnin = 3333;
rhos = zeros(1,num_of_samples-burnin);
s2s = zeros(1,num_of_samples-burnin); 
s_diagonal_terms = zeros(1,num_of_samples-burnin); 

for i=1:(num_of_samples-burnin)
    [A, Lambda] = uncollapseParameters(samples(:,i+burnin));
%     [d,~] = size(Lambda);
    U = (eye(n_timepoints) - A)/(eye(n_timepoints) + A);
    sigma = U*(Lambda.^2)*(U');
%     rho_av = 0;
%     for j = 0:58
%         rho_av = rho_av + sigma(1+j,2+j)/sigma(1+j,1+j);
%     end
%     rhos(i) = rho_av/59;
    rhos(i) = sigma(1,2)/sigma(1,1);
    s2s(i) = sigma(1,1)*(1 - (sigma(1,2)/sigma(1,1)).^2);
    s_diagonal_terms(i) = sigma(1,1);
end


upper95_rho = quantile(rhos,0.975);
lower5_rho = quantile(rhos,0.025);

lowers = mean(rhos >= rho);
uppers = mean(rhos <= rho);

iteration = str2num(env_id);
double_side_CI_coverage = (lower5_rho<=rho)&(rho<=upper95_rho);
table1 = [double_side_CI_coverage, lowers, uppers, iteration];
my_table = table(table1);
my_file_name = strcat("fixed_rho_outputs/arraynumber_", env_id, ".csv");
writetable(my_table, my_file_name);
    

upper95_sd = quantile(sqrt(s2s),0.975);
lower5_sd = quantile(sqrt(s2s),0.025);

lowers = mean(sqrt(s2s) >= s_true);
uppers = mean(sqrt(s2s) <= s_true);

iteration = str2num(env_id);
double_side_CI_coverage = (lower5_sd<=s_true)&(s_true<=upper95_sd);
table1 = [double_side_CI_coverage, lowers, uppers, iteration];
my_table = table(table1);
my_file_name = strcat("fixed_sd_outputs/arraynumber_", env_id, ".csv");
writetable(my_table, my_file_name);

upper975_sigma = quantile(s_diagonal_terms,0.975);
lower25_sigma = quantile(s_diagonal_terms,0.025);

diagonal_true = s_true.^2/(1-rho.^2);
lowers = mean(s_diagonal_terms >= diagonal_true)
uppers = mean(s_diagonal_terms <= diagonal_true);

iteration = str2num(env_id);
double_side_CI_coverage = (lower25_sigma<=diagonal_true)&(diagonal_true<=upper975_sigma);
table1 = [double_side_CI_coverage, lowers, uppers, iteration];
my_table = table(table1);
my_file_name = strcat("fixed_sigma_outputs/CIs_arraynumber_", env_id, ".csv");
writetable(my_table, my_file_name);


