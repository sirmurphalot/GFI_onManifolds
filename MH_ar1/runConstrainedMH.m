function [samples, accepts] = runConstrainedMH(num_iters, burn_in, ...
    proposal_scale, initial_guess, ll_function, ...
    constraintFunc, dConstraintFunc)
    [parameter_dim,~] = size(initial_guess);
    samples = zeros(parameter_dim, num_iters+1);
    accepts = zeros(1, num_iters);
    MH_ratios = zeros(1, 1);
    proposal_contribution = zeros(1, 1);
    likelihood_contribution = zeros(1, 1);
    likelihoods = zeros(1, 1);
    forward_probs = zeros(1, 1);
    backwards_probs = zeros(1, 1);
    flags = zeros(1,num_iters);
    samples(:,1)=initial_guess;
    
    dim=length(initial_guess);
    xx=dim;
	dd=(1/2)*(-1+sqrt(1+8*xx));
    
    total_iters = num_iters + burn_in;
    for iter=1:total_iters
        if mod(iter, 4)==0
            formatSpec = 'Acceptance rate (so far): %d.\n \n';
            str = sprintf(formatSpec,sum(accepts)/iter);
            disp(str); %pause(0.1);
        end
        % Calculate a proposal distribution by moving along the tangent
        % plane.
        disp("check the proposal value");
        n_timepoints = 10;
        > [tempA,tempL] = uncollapseParameters(samples(:,1));
        U = (eye(n_timepoints) - tempA)/(eye(n_timepoints) + tempA);
        sigma = U*(tempL.^2)*(U');
        [r,s] = AandLambdaToRhoandSigma(samples(:,1))
        ar1_constraint(y)
        [proposal, logdens_forward, ~, flag]=findProposal(samples(:,iter),...
            proposal_scale, constraintFunc, dConstraintFunc);
        % If this simulation moves off the manifold, automatically reject.
        
        if flag == 0
            fprintf(1,'\b\b\b%3.0i',iter);
            samples(:,iter+1)=samples(:,iter);
            continue
        end

	% Calculate the probability of moving back, on the tangent plane at
	% the proposal value.
	[logdens_reverse, flag] = findReverseProposalProb(samples(:,iter), ...
	    proposal, proposal_scale, constraintFunc, dConstraintFunc);
	% If the solver is unable to get back, reject automatically (for
	% reversability).
	if flag == 0
	    fprintf(1,'\b\b\b%3.0i',iter);
	    samples(:,iter+1)=samples(:,iter);
	    continue
	end
	flags(iter) = 1;
	
	% Put in a catch that makes sure that my proposal does not move off
	% of the hypercube.
	a_elements = proposal(1:(xx-dd));

	% Print out rho (debug) and make sure it is a valid rho.
	disp("sampled rho is");
	[r,~] = AandLambdaToRhoandSigma(proposal);
        [curr_r,~] = AandLambdaToRhoandSigma(samples(:,iter));
	if r>1 || r<-1
           samples(:,iter+1)=samples(:,iter);
           continue
        end
        
        if r<1
            disp("proposal rho is");
            disp(r);
            disp("current rho is");
            disp(curr_r);
        end
        
        [r,~] = AandLambdaToRhoandSigma(samples(:,iter));
        
        MH_ratio = ll_function(proposal) + logdens_reverse - ll_function(samples(:,iter)) - logdens_forward;
       
        if isnan(MH_ratio)
            samples(:,iter+1)=samples(:,iter);
            disp("MH ratio is not a number!");
            continue
        end
 
        % Accept or reject.
        
        disp("check the proposal value");
        
        %[~,L] = uncollapseParameters(samples(:,iter))
        

        unif_value = log(rand);
        if unif_value <= MH_ratio
            samples(:,iter+1)=proposal;
            accepts(:,iter)=1;
            MH_ratios(end+1) = MH_ratio;
            likelihoods(end+1) = ll_function(proposal);
            
            forward_probs(end+1) = logdens_forward;
            backwards_probs(end+1) = logdens_reverse;
            
            proposal_contribution(end+1) = logdens_reverse - logdens_forward;
            likelihood_contribution(end+1) = ll_function(proposal) - ll_function(samples(:,iter));
        else
            samples(:,iter+1)=samples(:,iter);
            MH_ratios(end+1) = MH_ratio;
            likelihoods(end+1) = ll_function(samples(:,iter));
            
            forward_probs(end+1) = logdens_forward;
            backwards_probs(end+1) = logdens_reverse;
            
            proposal_contribution(end+1) = logdens_reverse - logdens_forward;
            likelihood_contribution(end+1) = ll_function(proposal) - ll_function(samples(:,iter));
        end
        
        disp(proposal_contribution(end));
        

	if mod(iter,100)==0
% 		temp_samples = samples(:,1:(iter+1));
% 		save("data_logs/sample_data_running.mat",'temp_samples');
%         	temp_MH = MH_ratios;
% 		save("data_logs/MH_ratios_running.mat",'temp_MH');
%         	temp_likelihood = likelihoods;
% 		save("data_logs/likelihoods_running.mat",'temp_likelihood');
%         
%         	temp_forward = forward_probs;
% 		save("data_logs/forward_probs_running.mat",'temp_forward');
%         	temp_backward = backwards_probs;
% 		save("data_logs/backward_probs_running.mat",'temp_backward');
%         
%         	temp_proposal = proposal_contribution;
% 		save("data_logs/proposal_contributions_running.mat",'temp_proposal');
%         	temp_likelihood_contribution = likelihood_contribution;
% 		save("data_logs/likelihood_contributions_running.mat",'temp_likelihood_contribution');
	end        
        fprintf(1,'\b\b\b%3.0i',iter);
    end
    samples = samples(:,(burn_in+1):(num_iters+1));
    accepts = accepts(:,(burn_in):num_iters);
end
