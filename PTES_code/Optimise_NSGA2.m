clear all
clc
global V M xl xu etac etam p pop_size pm

%% Description

% 1. This is the main program of NSGA II. It requires only one input, which is test problem
%    index, 'p'. NSGA II code is tested and verified for 14 test problems.
% 2. This code defines population size in 'pop_size', number of design
%    variables in 'V', number of runs in 'no_runs', maximum number of 
%    generations in 'gen_max', current generation in 'gen_count' and number of objectives
%    in 'M'.
% 3. 'xl' and 'xu' are the lower and upper bounds of the design variables.
% 4. Final optimal Pareto soutions are in the variable 'pareto_rank1', with design
%    variables in the coumns (1:V), objectives in the columns (V+1 to V+M),
%    constraint violation in the column (V+M+1), Rank in (V+M+2), Distance in (V+M+3).
%% code starts
M=2;                     %Number of objectives
pop_size=5;             % Population size
no_runs=1;               % Number of runs
gen_max=5;              % MAx number of generations - stopping criteria
 
V=3;                     %Nr of variables
xl=[0   0.8  0.9];      % lower bound vector
xu=[0   1.6  0.97];     % upper bound vector
         

etac = 20;                  % distribution index for crossover
etam = 20;                  % distribution index for mutation / mutation constant

pm=1/(V);                     % Mutation Probability

fname='PTES_opt';
Q=[];
for run = 1:no_runs
    
%% Initial population 
xl_temp=repmat(xl, pop_size,1);
xu_temp=repmat(xu, pop_size,1);
x = xl_temp+((xu_temp-xl_temp).*rand(pop_size,V));

%% Evaluate objective function
for i =1:pop_size
    try
    [ff(i,:) err(i,:)] =feval(fname, x(i,:));          % Objective function evaulation 
    continue;
    end
end
error_norm=normalisation(err);                      % Normalisation of the constraint violation
population_init=[x ff error_norm];
[population front]=NDS_CD_cons(population_init); % Non domination Sorting on initial population
    
%% Generation Starts
for gen_count=1:gen_max
% selection (Parent Pt of 'N' pop size)
parent_selected=tour_selection(population);                     % 10 Tournament selection
%% Reproduction (Offspring Qt of 'N' pop size)
child_offspring  = genetic_operator(parent_selected(:,1:V));    % SBX crossover and polynomial mutation
if size(child_offspring,1) ~= pop_size %Issue arises when rows with same values are discarded by the genetic operator. Usually happens for only 1 row.
    [C,ia,ic] = unique(parent_selected,'rows') ;
    repeatedrow=numel(ic);
    child_offspring(pop_size,1:V)=parent_selected(repeatedrow-1,1:V);
    
end

for ii = 1:pop_size
    try
        [fff(ii,:) err(ii,:)]=feval(fname, child_offspring(ii,:));      % objective function evaluation for offspring
    catch
    continue;
    end
end

error_norm=normalisation(err);                                  
child_offspring=[child_offspring fff error_norm];

%% INtermediate population (Rt= Pt U Qt of 2N size)
population_inter=[population(:,1:V+M+1) ; child_offspring(:,1:V+M+1)];
[population_inter_sorted front]=NDS_CD_cons(population_inter);              % Non domination Sorting on offspring
%% Replacement - N
new_pop=replacement(population_inter_sorted, front);
population=new_pop;
end
new_pop=sortrows(new_pop,V+1);
paretoset(run).trial=new_pop(:,1:V+M+1);
Q = [Q; paretoset(run).trial];                      % Combining Pareto solutions obtained in each run
end

%% Result and Pareto plot
if run==1
for iii =1:pop_size %This loop is a repetition just to find the capital cost. Needs replacement
    try
    [final(iii,:) err(iii,:) extra(iii,:)] =feval(fname, new_pop(iii,:));          % Objective function evaulation 
    continue;
    end
end
format long
Final_results=[new_pop(:,1:V+M) extra]
plot(new_pop(:,V+1),new_pop(:,V+2),'*')
else                                        
[pareto_filter front]=NDS_CD_cons(Q);               % Applying non domination sorting on the combined Pareto solution set
rank1_index=find(pareto_filter(:,V+M+2)==1);        % Filtering the best solutions of rank 1 Pareto
pareto_rank1=pareto_filter(rank1_index,1:V+M);
plot(pareto_rank1(:,V+1),pareto_rank1(:,V+2),'*')   % Final Pareto plot
end
xlabel('1-Efficiency')
ylabel('LCOS')
title(' PTES')