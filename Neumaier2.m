function f=Neumaier2(x)
% Neumaier 2 Problem (NF2) n=4
% The number of local minima is unknown
% LimInf=[0 0 0 0]; LimSup=[n n n n];
% Fojmin=0
% x*=[1 2 2 3]
n=length(x); b=[8 18 44 114]; 
S1=zeros(1,n); S3=S1;
for k=1:n
    for i=1:n
        S1(i)=x(i)^k;
    end
    S2=sum(S1);
    S3(k)=(b(k)-S2)^2;
end
f=sum(S3);
end