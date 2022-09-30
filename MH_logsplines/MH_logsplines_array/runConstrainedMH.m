function [samples, accepts] = runConstrainedMH(num_iters, burn_in, ...
    proposal_scale, initial_guess, ll_function, ...
    constraintFunc, dConstraintFunc)
% Author: Alexander Murph
    [parameter_dim,~] = size(initial_guess);
    samples = zeros(parameter_dim, num_iters+1);
    accepts = zeros(1, num_iters);
    samples(:,1)=initial_guess;
    
    total_iters = num_iters + burn_in;
    
    for iter=1:total_iters
        if mod( iter , 250 ) == 0
            formatSpec = 'Finished %d iterations.';
            str = sprintf(formatSpec,iter);
            disp(str);
            formatSpec = 'Calculated %d accepts.';
            str = sprintf(formatSpec,sum(accepts));
            disp(str);
            formatSpec = 'Current value: %d %d %d %d % %d %d %d %d %d.';
            str = sprintf(formatSpec,samples(:,iter));
            disp(str);
        end
        [proposal, prob_forward, Tx, flag]=findProposal(samples(:,iter),...
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
        if iter > burn_in
%              disp("hey");
        end
        MH_ratio = ll_function(proposal) + log(prob_reverse) - ll_function(samples(:,iter)) - log(prob_forward);
        unif_value = log(rand);
        if isinf(MH_ratio)
            disp("hey");
        end
        
        if mod( iter , 100 ) == 0
            formatSpec = 'The current MH ratio is: %d';
            str = sprintf(formatSpec,MH_ratio);
            disp(str);
        end
        
        if unif_value <= MH_ratio
            samples(:,iter+1)=proposal;
            accepts(:,iter)=1;
            continue
        end
        samples(:,iter+1)=samples(:,iter);
    end
%     samples = samples(:,(burn_in+1):(num_iters+1));
%     accepts = accepts(:,(burn_in):num_iters);
end