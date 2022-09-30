env_id = getenv(char("SLURM_ARRAY_TASK_ID"))
%env_id = '2'
rng(str2num(env_id))

n_timepoints = 10;
n_observations = 1;
data_matrix = zeros(n_observations, n_timepoints);
s_true = 1
rho = 0.7
sigma_true = s_true.^2/(1-rho.^2)


first_value = normrnd(0,s_true/sqrt(1-rho.^2),1);
data_matrix(1) = first_value;
for i = 1:n_observations
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

disp("getting the conditional version")
rho_cond_mle = estimates(2,1)
sd2_cond_mle = 0
for i = 2:n_timepoints
    sd2_cond_mle = sd2_cond_mle + ((data_matrix(1,i) - estimates(1,1) - ...
        rho_cond_mle*data_matrix(1,i-1)).^2)/(n_timepoints-1);
end
sd2_cond_mle 
my_ll_func = @(rs)(-ll(rs(1),rs(2),data_matrix(1,:)));

rs_guess = [rho_cond_mle sqrt(sd2_cond_mle)];
rs_max = fminsearch(my_ll_func, rs_guess);

disp("getting the numeric maximization version")
rho_cond_mle = rs_max(1)
sd2_cond_mle = rs_max(2)
rho_cond_mle = 0.7;
sd2_cond_mle=1;
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
disp("total number of matrices is")
disp(2^d-1)

% We do a preliminary run to find a single skew-symmetric matrix that works
for i = 0:(2^d-1)
    if mod( i , 5000 ) == 0
        disp(i);
    end
    temp_signs = de2bi(i,d);
    temp_signs = 1 - 2*temp_signs(1:d);
    % By O'Dorney 2014, signature matrices with an odd number of -1's will
    % not have a corresponding Cayley Transform.
    if(1==original_sign*prod(temp_signs) || ( mod(sum(temp_signs == -1),2)==1 ))
        tempUhat = U*diag(temp_signs);
        
        % We double check that the I+U matrix is non-singular (it must be).
        is_nonsingular = eye(d) + tempUhat;
        if(det(is_nonsingular)==0)
            continue
        else 
            skew_symmetric_matrix = (eye(d) - tempUhat)/(eye(d) + tempUhat);
	    break
        end
    end
end

U = tempUhat;
original_sign = sign(det(U));

% O'Dorney 2014 says that the signature matrix that we want is the one that
% maximizes H_D(S).
 U = tempUhat;
 current_maximum = -1;
 for i = 0:(2^d-1)
     if mod( i , 5000 ) == 0
         disp(i);
     end
     temp_signs = de2bi(i,d);
     temp_signs = 1-2*temp_signs(1:d);
     % By O'Dorney 2014, signature matrices with an odd number of -1's will
     % not have a corresponding Cayley Transform.
     if(1==original_sign*prod(temp_signs) )
         tempUhat = U*diag(temp_signs);
         
         % We double check that the I+U matrix is non-singular (it must be).
         is_nonsingular = eye(d) + tempUhat;
         if(det(is_nonsingular)==0)
             continue
         end
         
         new_hidden_value = find_hidden_entry(temp_signs, skew_symmetric_matrix);
          if(current_maximum > new_hidden_value)% %              disp("found new max");
 %              disp("found new max");
              continue
          else 
              current_maximum = new_hidden_value;
          end
	
%          E = eig(tempUhat);
%          if ~((sum(abs(real(E)+1)<1e-4)==0) )
%              continue
%          end         
         
         tempA = (eye(d) - tempUhat)/(eye(d) + tempUhat);
         temp_vector = collapseParameters(tempA,Lambda);
 
         xx = d*(d+1)/2;
         a_elements = temp_vector(1:(xx-d));
         if((1 == prod(a_elements>-1 & a_elements<1)) )
             disp("got here");
             [actualA,~] = uncollapseParameters(temp_vector);
             actualU = tempUhat;
             %initial_guess = temp_vector;
             found_starting_point = true;
             break
         end
         
 %         if ~(sum(abs(real(E)+1)<1e-4)==0)
 %             count = count + 1;
 %             continue;
 %         else
 %             disp("found one!")
 %             U = tempUhat;
 %             break
 %         end
     end
 end

% if(found_starting_point)
%     disp("OH MY GOODNESS IT FOUND ONE YAY!");
%     my_file_name_data = strcat("working_A_mats/temp_vector_", env_id, ".mat");
% 
% 
% 	save(my_file_name_data,'temp_vector');
% end


% A = [];
% b = [];
% Aeq = [];
% beq = [];
% lb = zeros([1,d]);
% ub = ones([1,d]);
% f = @(a)(minimizeHiddenFunction(a, skew_symmetric_matrix));
% fun = @(a)(-find_hidden_entry(a, skew_symmetric_matrix));
% intcon = 1:d;
% initial = zeros([1,d]);
%i = floor(2^d/4+3^(floor(d/2)));
%initial = de2bi(i,d);

%options = optimoptions('ga','ConstraintTolerance',1e-30, 'FunctionTolerance', 1e-30, 'FunctionTolerance', 1e-20,...
%                        'MaxGenerations', 1e5, 'MaxStallGenerations', 1e5);

%my_const = @(a)(optim_constraint(a));
                     
%x = ga(fun,d,A,b,Aeq,beq,lb,ub,my_const,intcon,options);


%disp("hey");
%temp_signs = x;
%temp_signs = 1 - 2*temp_signs(1:d);
%tempUhat = tempUhat*diag(temp_signs);

% new_hidden_value = find_hidden_entry(temp_signs, skew_symmetric_matrix);


disp("creating A matrix");
A = (eye(n_timepoints)-tempUhat)/(eye(n_timepoints)+tempUhat);

% Note that this will move it off the manifold, but im willing to try.
%A(find(A<-1)) = -1;
%A(find(A>1)) = 1;
% A = A_guess;
disp("created A matrix");

initial_guess = collapseParameters(A,Lambda);
% initial_guess = [A';diag(Lambda)];
% [A,Lambda] = uncollapseParameters(initial_guess);

% Notice that 
starting_covariance(2,1)/starting_covariance(1,1)

% Yet
disp("trying to make U matrix again");
newU = (eye(n_timepoints)-A)/(eye(n_timepoints)+A);
newCov = newU * (Lambda.^2) * (newU');
newCov(2,1)/newCov(1,1)
disp("made U matrix.");

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
num_iters = 500000;
burn_in = 1;
% paper used 2e-3 as of Sep 7th
proposal_scale = 0.07

ll_function = @(x) ll_density(x, data_matrix);
constraintFunc = @(x) (ar1_constraint(x));
dConstraintFunc = @(x) (dar1_constraint(x));
disp("Entering the MCMC Algorithm.")
warning('off');
disp("starting the MCMC sampling.")
[samples, accepts] = runConstrainedMH(num_iters, burn_in, ...
    proposal_scale, initial_guess, ll_function, ...
    constraintFunc, dConstraintFunc);
mean(accepts)


my_file_name_data = strcat("data_logs/sample_data_", env_id, ".mat");


save(my_file_name_data,'samples');
%samples = load(my_file_name_data);
%samples = samples.samples;
%% Use the samples to calculate values for rho and sigma2:


[~,num_of_samples] = size(samples);
burnin = 150000;
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


upper95_rho = quantile(rhos,0.975)
lower5_rho = quantile(rhos,0.025)


%lowers = mean(rhos >= rho_cond_mle)
%uppers = mean(rhos <= rho_cond_mle)
lowers = mean(rhos >= rho)
uppers = mean(rhos <= rho)

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

%lowers = mean(sqrt(s2s) >= sqrt(sd2_cond_mle));
%uppers = mean(sqrt(s2s) <= sqrt(sd2_cond_mle));

iteration = str2num(env_id);
double_side_CI_coverage = (lower5_sd<=s_true)&(s_true<=upper95_sd);
table1 = [double_side_CI_coverage, lowers, uppers, iteration];
my_table = table(table1);
my_file_name = strcat("fixed_sd_outputs/arraynumber_", env_id, ".csv");
writetable(my_table, my_file_name);

upper975_sigma = quantile(s_diagonal_terms,0.975)
lower25_sigma = quantile(s_diagonal_terms,0.025)

diagonal_true = s_true.^2/(1-rho.^2);
lowers = mean(s_diagonal_terms >= diagonal_true)
uppers = mean(s_diagonal_terms <= diagonal_true);

iteration = str2num(env_id);
double_side_CI_coverage = (lower25_sigma<=diagonal_true)&(diagonal_true<=upper975_sigma);
table1 = [double_side_CI_coverage, lowers, uppers, iteration];
my_table = table(table1);
my_file_name = strcat("fixed_sigma_outputs/CIs_arraynumber_", env_id, ".csv");
writetable(my_table, my_file_name);



