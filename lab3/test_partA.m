ObjectiveFunction = @mysum;
nvars = 20;    % Number of variables
LB = zeros(1,20);   % Lower bound
UB = ones(1,20);  % Upper bound
%ConstraintFunction = @simple_constraint;
[x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB)
