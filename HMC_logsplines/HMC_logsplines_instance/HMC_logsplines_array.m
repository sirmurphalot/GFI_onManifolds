env_id = getenv(char("SLURM_ARRAY_TASK_ID"))
%env_id = "1"
rng(str2num(env_id))

n=500;

lower = 0;
peak = 0.2;
upper = 1;
pd = makedist('Triangular','A',lower,'B',peak,'C',upper);
y = random(pd,1,n);

knots = [-0.0250:0.1:1.1];
knots(end-1) = 1;
knots(2) = 0;
constraintFunc = @(x) logspline_constraint(knots, x);
options = optimoptions('fsolve','Display','none','TolFun', 1e-10);

flag=11;
solve_this = @(a) (constraintFunc(a));
[f,xi] = ksdensity(y,knots);
initial = [log(f(2)), log(f(3)), log(f(4)), log(f(5)), log(f(6)), log(f(7)),...
            log(f(8)), log(f(9)), log(f(10)), log(f(11))];
q0 = fsolve(solve_this, initial, options);


nllFunc = @(x) negLog_fiducial_likelihood(knots,x,y);
conFunc = constraintFunc;

% Set parameters (that I don't quite understand)
nRuns = 1+10;
Nscale = 100000;
runTime = 15; % Time per run in seconds
nTrials = 100;
doPrint = false;

samplers = struct([]);

h = 2e-5;

M = 1;
L = 5;
hmcOpts = struct('intMethod',0,'doPrint',doPrint,'printMod',10);
samplers(end+1).name = sprintf('CHMC (L = %d)',L);
samplers(end).desc = sprintf('M = %g, h = %g',M,h);
samplers(end).func = @(q,N) constrainedHMC(q,nllFunc,conFunc,M,N,L,h,hmcOpts);
samplers(end).N = ceil(Nscale./L);
samplers(end).style = 'b-';

samplerQs = cell(length(samplers),1+nTrials);
samplerStats = cell(length(samplers),1+nTrials);
samplerNLLs = cell(length(samplers),1+nTrials);
samplerTs = zeros(length(samplers),1+nTrials);
for i = 1:numel(samplers)
    fprintf('%s...',samplers(i).name);
    [q1,qs,samplerStats{i,1}] = samplers(i).func(q0,samplers(i).N);
    samplerQs{i,1} = [q0' qs q1];
    samplerTs(i,1) = samplerStats{i,1}.t/samplerStats{i,1}.N;
%     samplerNLLs{i,1} = nllFunc(samplerQs{i,1});
    
    currN = ceil(runTime./samplerTs(i,1));
    fprintf('%d samples: ',currN);
    fprintf('\n');
end
mean(samplerStats{1,1}.accepted)
Qs = samplerQs{1,1};
% samples = Qs.Qs;

[num_rows, num_cols] = size(samples);
num_knots = length(knots) - 2;
dens_values = zeros(num_knots,num_cols);
upperCIs = zeros(1,num_knots);
lowerCIs = zeros(1,num_knots);
CI_coverages = zeros(1,num_knots);


% This is where I am going to need to process the simulation values.  
% First I'll get the upper and lower bounds on the logspline density values
% AT the knots.
for i=1:num_knots
    for j=1:num_cols
        dens_values(i,j) = exp(logspline_density(knots(i+1),knots,samples(:,j)));
    end
end

for i=1:num_knots
%     uppers = quantile(dens_values(i,:), 0.975);
%     lowers = quantile(dens_values(i,:), 0.025);
    density_of_triangle_distn = pdf(pd,knots(i+1));
    lowerCIs(i) = mean(density_of_triangle_distn >= dens_values(i,:));
    upperCIs(i) = mean(density_of_triangle_distn <= dens_values(i,:));
    
    uppers = quantile(dens_values(i,:), 0.975);
    lowers = quantile(dens_values(i,:), 0.025);
    
    CI_coverages(i) = (density_of_triangle_distn<=uppers)&(density_of_triangle_distn>=lowers);
end

iteration = str2num(env_id);
table1 = [lowerCIs, upperCIs, CI_coverages];
T = array2table(table1);
T.Properties.VariableNames(1:(num_knots*3)) = {'lowerCI_knot2', 'lowerCI_knot3', 'lowerCI_knot4',...
                                                 'lowerCI_knot5', 'lowerCI_knot6', 'lowerCI_knot7', 'lowerCI_knot8',...
                                                 'lowerCI_knot9', 'lowerCI_knot10', 'lowerCI_knot11',...
                                                 'upperCI_knot2', 'upperCI_knot3', 'upperCI_knot4',...
                                                 'upperCI_knot5', 'upperCI_knot6', 'upperCI_knot7', 'upperCI_knot8',...
                                                 'upperCI_knot9', 'upperCI_knot10', 'upperCI_knot11',...
                                                 'CI_coverage2', 'CI_coverage3', 'CI_coverage4', 'CI_coverage5','CI_coverage6',...
                                                 'CI_coverage7', 'CI_coverage8', 'CI_coverage9', 'CI_coverage10','CI_coverage11'};
my_file_name = strcat("CI_outputs/arraynumber_", env_id, ".csv");
writetable(T,my_file_name);


