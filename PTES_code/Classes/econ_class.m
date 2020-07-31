classdef econ_class
   properties
       % Costs
       cost_mode
       COST     % Total capital cost, $
       cost     % Specific capital cost - e.g. Cost / power
       
       sd      % Standard deviation (for uncertainty analysis) (fraction of mean)
       ul      % upper cost limit (fraction of mean)
       ll      % lower cost limit (fraction of mean)
       
   end
   
   methods
       function obj = econ_class(cost_mode, sd, ul, ll)
            
           obj.cost_mode = cost_mode ;
           obj.sd = sd ;
           obj.ul = ul ;
           obj.ll = ll ;
            
       end
       
       % Calculate distribution of costs
       function [cost_array] = cost_sens(obj, N)
           
           % N is the number of points to evaluate
           MN = obj.COST ;
           SD = obj.sd * MN ; % Standard deviation
           LL = obj.ll * MN * ones(N,1) ; % Lower bound
           UL = obj.ul * MN * ones(N,1) ; % Upper bound
           
           cost_array = trandn((LL-MN)./SD,(UL-MN)./SD) ; % Truncated normal distribution array
           cost_array = MN + SD .* cost_array ;            % Transform
                      
       end
   end
end