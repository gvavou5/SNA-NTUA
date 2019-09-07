%save('optimproblemdouble.mat','optimproblem');
%save('optimresultsdouble.mat','optimresults');
%save('optionsdouble.mat','options');
%save('optimproblem.mat','optimproblem');
%save('optimresults.mat','optimresults');
%save('options.mat','options');


% Load these if you want to use Bit String 
load('optimproblem.mat');
load('optimresults.mat');
load('options.mat');

% Load these if you want to use Double
%load('optimproblemdouble.mat');
%load('optimresultsdouble.mat');
%load('optionsdouble.mat');

indx1 = 0;
indx2 = 0;
indx3 = 0;

% i load to workspace the parameteres via optimization app and the change
% with the following for loops the requested fields. If i want to do it
% manually i drop the comments below.

%ObjectiveFunction = @mysum;
%nvars = 20;         % Number of variables
%LB = zeros(1,20);   % Lower bound
%UB = ones(1,20);    % Upper bound

for i = 10:10:100
    indx2=indx2+1;
    for j = 0.3:0.1:0.9
        indx3=indx3+1;
        for y = 0.01:0.01:0.2
            indx1=indx1+1;
            %optimproblem.solver = 'ga';
            %optimproblem.options = optimoptions('ga');
            %optimproblem.options.PopulationType = 'bitstring';
            optimproblem.options.PopulationSize = i;
            optimproblem.options.CrossoverFraction = j;
            optimproblem.options.MutationFcn{2} = y;
            %[x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB);
            [x,fval] = ga(optimproblem);
            A(indx1,indx3,indx2)=sum(x);
        end
        indx1=0;
    end
    indx3=0;
end

%save('Part1.mat','A');

mkdir Part1A;
for n = 10:10:100
    figure(n/10);
    bar3(A(:,:,n/10));
    X=['Plot of the sum of x1,x2...,xn for Popoulation = ' num2str(n)];
    title(X);
    cd Part1A; print -djpeg; cd ..
end
close all;
