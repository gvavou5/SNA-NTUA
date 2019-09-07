function value = ego_net(A,n)
%EGO_NET Summary of this function goes here
%   Detailed explanation goes here

%B=A*A;
%check=triu(1-A,1);
%sum=0.0;
%for  i=1:n 
%    for j=1:n
%        if check(i,j) ==1
%            sum=sum+1/B(i,j);
%        end

    %end
    %end    
%value =sum;
t=brandesBetwCentr(A);
value=t(1);