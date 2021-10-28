function [samples, accepts] = runConstrainedMH(num_iters, burn_in, ...
    proposal_scale, initial_guess, ll_function, ...
    constraintFunc, dConstraintFunc)
    [parameter_dim,~] = size(initial_guess);
    samples = zeros(parameter_dim, num_iters+1);
    accepts = zeros(1, num_iters);
    flags = zeros(1,num_iters);
    samples(:,1)=initial_guess;
    
    total_iters = num_iters + burn_in;
    
    for iter=1:total_iters
        [proposal, prob_forward, ~, flag]=findProposal(samples(:,iter),...
            proposal_scale, constraintFunc, dConstraintFunc);
        if flag == 0
            samples(:,iter+1)=samples(:,iter);
            continue
        end
        
        [prob_reverse, flag] = findReverseProposalProb(samples(:,iter), ...
            proposal, proposal_scale, constraintFunc, dConstraintFunc);
        if flag == 0
            samples(:,iter+1)=samples(:,iter);
            continue
        end
        flags(iter) = 1;
        if iter > burn_in
%             disp("hey");
        end
        MH_ratio = ll_function(proposal) + log(prob_reverse) - ...
            ll_function(samples(:,iter)) - log(prob_forward);
        formatSpec = 'MH ratio is: %d.';
        str = sprintf(formatSpec,MH_ratio);
        disp(str);
        if isinf(MH_ratio)
            disp("hey");
        end
        
        unif_value = log(rand);
        
        if unif_value <= MH_ratio
            samples(:,iter+1)=proposal;
            accepts(:,iter)=1;
            continue
        end
        samples(:,iter+1)=samples(:,iter);
    end
    samples = samples(:,(burn_in+1):(num_iters+1));
    accepts = accepts(:,(burn_in):num_iters);
end