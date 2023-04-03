env_id = getenv(char("SLURM_ARRAY_TASK_ID"))
rng(str2num(env_id))

my_file_name_data = strcat("data_logs/sample_data_", env_id, ".mat");


samples = load(my_file_name_data);
samples = table(samples.samples);
%% Use the samples to calculate values for rho and sigma2:
my_file_name = strcat("full_data_csvs/data_from_arraynumber_", env_id, ".csv")
writetable(samples, my_file_name);


