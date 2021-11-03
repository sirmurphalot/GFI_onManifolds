function [samples, accepts] = runConstrainedMH(num_iters, burn_in, ...
    proposal_scale, initial_guess, ll_function, ...
    constraintFunc, dConstraintFunc)
    [parameter_dim,~] = size(initial_guess);
    [matrix_dim,~] = size(uncollapseParameters(initial_guess));
    samples = zeros(parameter_dim, num_iters+1);
    accepts = zeros(1, num_iters);
    flags = zeros(1,num_iters);
    samples(:,1)=initial_guess;
    
    total_iters = num_iters + burn_in;
    for iter=1:total_iters
        % Calculate a proposal distribution by moving along the tangent
        % plane.
        [proposal, logdens_forward, ~, flag]=findProposal(samples(:,iter),...
            proposal_scale, constraintFunc, dConstraintFunc);
        % If this simulation moves off the manifold, automatically reject.
        if flag == 0
            samples(:,iter+1)=samples(:,iter);
            continue
        end
%         % I'm also going to constrain the values of A so that they always
%         % remain on a hypercube between -1 and 1 (see Murph et al 2021).
%         temp_vector = proposal;
%         if(~(1 == prod(temp_vector((matrix_dim+1):end)>-1 & temp_vector((matrix_dim+1):end)<1)) )
%             disp("Proposal moved off of the hypercube!");
%             samples(:,iter+1)=samples(:,iter);
%             continue
%         end
        
        % Calculate the probability of moving back, on the tangent plane at
        % the proposal value.
        [logdens_reverse, flag] = findReverseProposalProb(samples(:,iter), ...
            proposal, proposal_scale, constraintFunc, dConstraintFunc);
        % If the solver is unable to get back, reject automatically (for
        % reversability).
        if flag == 0
            samples(:,iter+1)=samples(:,iter);
            continue
        end
        flags(iter) = 1;
        
        % Calculate the Metropolis ratio.
        MH_ratio = ll_function(proposal) + logdens_reverse - ...
            ll_function(samples(:,iter)) - logdens_forward;
        formatSpec = 'MH ratio is: %d.';
        str = sprintf(formatSpec,MH_ratio);
        disp(str);
        formatSpec = '----> Portion from the likelihood: %d.';
        str = sprintf(formatSpec,ll_function(proposal)-ll_function(samples(:,iter)));
        disp(str);
        formatSpec = '----> Portion from the proposals: %d.';
        str = sprintf(formatSpec,logdens_reverse-logdens_forward);
        disp(str);
        formatSpec = '--------> (neg) Log density forward: %d.';
        str = sprintf(formatSpec,-logdens_forward);
        disp(str);
        formatSpec = '--------> Log density backward: %d.';
        str = sprintf(formatSpec,logdens_reverse);
        disp(str);
        
        % Accept or reject.
        unif_value = log(rand);
        if unif_value <= MH_ratio
            samples(:,iter+1)=proposal;
            accepts(:,iter)=1;
        else
            samples(:,iter+1)=samples(:,iter);
        end
        
        formatSpec = 'Acceptance rate (so far): %d.';
        str = sprintf(formatSpec,sum(accepts)/iter);
        disp(str);
        formatSpec = "Iteration number is: %i.";
        str = sprintf(formatSpec,iter);
        disp(str);
        disp("-----------");
        
    end
    samples = samples(:,(burn_in+1):(num_iters+1));
    accepts = accepts(:,(burn_in):num_iters);
end