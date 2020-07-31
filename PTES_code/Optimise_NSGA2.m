clear all
clc
global V M xl xu etac etam pop_size pm

%% Description

% 1. This code defines population size in 'pop_size', number of design
%    variables in 'V', number of runs in 'no_runs', maximum number of 
%    generations in 'gen_max', current generation in 'gen_count' and number of objectives
%    in 'M'.
% 2. 'xl' and 'xu' are the lower and upper bounds of the design variables.
% 3. Final optimal Pareto soutions are in the variable 'pareto_rank1', with design
%    variables in the coumns (1:V), objectives in the columns (V+1 to V+M),
%    constraint violation (if the input was provided) in the column (V+M+1), Rank in (V+M+2), Distance in (V+M+3).
%% code starts
M=2;                     %Number of objectives
pop_size=3;             % Population size
no_runs=1;               % Number of runs
gen_max=3;              % MAx number of generations - stopping criteria
 
V=3;                     %Nr of variables
xl=[523   0.8  0.8];      % lower bound vector
xu=[723   1.6  0.97];     % upper bound vector
         

etac = 20;                  % distribution index for crossover
etam = 20;                  % distribution index for mutation / mutation constant

pm=1/(V);                     % Mutation Probability

fname='PTES_opt';
Q=[];
x=zeros(pop_size,V);

Obj1='1 - Roundtrip efficiency';
Obj2='LCOS [USD/kWh]';
Obj3='Capital Cost [USD]';

Par1='Compressor Inlet Temperature';
Par2='Number of expansion stages';
Par3='Effectiveness';
for run = 1:no_runs
    
%% Initial population 
%Want to use previous data? change previousdata to 1
previousdata=0;

if previousdata==0
xl_temp=repmat(xl, pop_size,1);
xu_temp=repmat(xu, pop_size,1);
x = xl_temp+((xu_temp-xl_temp).*rand(pop_size,V));
end

%Getting previous data from the IntialPop.csv file in outputs folder
if previousdata==1
    A=csvread('./Data/POP');
    lengthcheck=size(x,1)-size(A,1);
    
    if lengthcheck==0
        x=A(:,1:V);
    end
    
    if lengthcheck>0
        x(1:size(A,1),1:V)=A(:,1:V);
        L1=size(A,1)+1; %if size(Nr of rows) of x is 40 and A is 30, then L1=31 and L2=40
        L2=size(x,1);
        for jj=L1:L2
            pop=xl+rand.*(xu-xl);
            x(jj,:)=pop;
        end
    end
        
     if lengthcheck<0
         x=A(1:pop_size,1:V);
     end      
end
    
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

if size(child_offspring,1) ~= pop_size   %Issue arises when rows with same values are discarded by the genetic operator when an error runs through.Usually happens for only 1 row.
    %child_offspring(pop_size,1:V)= child_offspring(randi(pop_size-1),1:V); %Randomly select a row in offspring to be added as the last row
    child_offspring(pop_size,1:V)= child_offspring(pop_size-1,1:V);        %Repeat the second last row as the last row
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
Final_results=[new_pop(:,1:V+M) extra];

if M==2
    figure(1)
    scatter(new_pop(:,V+1),new_pop(:,V+2),'*')
    hold on;
    grid on; xlabel(Obj1); ylabel(Obj2);
    title('PTES'); hold off;
    figure(2)
    scatter(Final_results(:,1),Final_results(:,4));
    grid on; xlabel(Par1); ylabel(Obj1);
    figure(3)
    scatter(Final_results(:,1),Final_results(:,5));
    grid on; xlabel(Par1); ylabel(Obj2);
    figure(4)
    scatter(Final_results(:,2),Final_results(:,4));
    grid on; xlabel(Par2); ylabel(Obj1);
    figure(5)
    scatter(Final_results(:,2),Final_results(:,5));
    grid on; xlabel(Par2); ylabel(Obj2);
    figure(6)
    scatter(Final_results(:,3),Final_results(:,4));
    grid on; xlabel(Par3); ylabel(Obj1);
    figure(7)
    scatter(Final_results(:,3),Final_results(:,5));
    grid on; xlabel(Par3); ylabel(Obj2);
end

if M==3
   figure(1)
   scatter3(new_pop(:,V+1),new_pop(:,V+2),new_pop(:,V+3),'ok'); 
   hold on;
   grid on; xlabel(Obj1); ylabel(Obj2); zlabel(Obj3);
   title(' PTES');
   hold off;
   figure(2)
   scatter(Final_results(:,1),Final_results(:,4));
   grid on; xlabel(Par1); ylabel(Obj1);
   figure(3)
   scatter(Final_results(:,1),Final_results(:,5));
   grid on;xlabel(Par1); ylabel(Obj2);
   figure(4)
   scatter(Final_results(:,1),Final_results(:,6));
   grid on;xlabel(Par1); ylabel(Obj3);
   figure(5)
   scatter(Final_results(:,2),Final_results(:,4));
   grid on; xlabel(Par2); ylabel(Obj1);
   figure(6)
   scatter(Final_results(:,2),Final_results(:,5));
   grid on; xlabel(Par2); ylabel(Obj2);
   figure(7)
   scatter(Final_results(:,2),Final_results(:,6));
   grid on; xlabel(Par2); ylabel(Obj3);
   figure(8)
   scatter(Final_results(:,3),Final_results(:,4));
   grid on; xlabel(Par3); ylabel(Obj1);
   figure(9)
   scatter(Final_results(:,3),Final_results(:,5));
   grid on; xlabel(Par3); ylabel(Obj2);
   figure(10)
   scatter(Final_results(:,3),Final_results(:,6));
   grid on; xlabel(Par3); ylabel(Obj3);
end


else                                        
[pareto_filter front]=NDS_CD_cons(Q);               % Applying non domination sorting on the combined Pareto solution set
rank1_index=find(pareto_filter(:,V+M+2)==1);        % Filtering the best solutions of rank 1 Pareto
pareto_rank1=pareto_filter(rank1_index,1:V+M);
plot(pareto_rank1(:,V+1),pareto_rank1(:,V+2),'*')   % Final Pareto plot
end

%csvwrite('./Data/IntialPop_solar',Final_results);