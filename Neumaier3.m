function f=Neumaier3(x)
% Neumaier 3 Problem (NF3) n=10
% The number of local minima is unknown
% LimInf=[-n^2 -n^2 ... -n^2]; LimSup=[n^2 n^2 ... n^2];
% Fojmin=-n(n+4)(n-1)/6
% xi*=i(n+1-i)
n=length(x); S1=zeros(1,n); S2=S1;
for i=1:n
    S1(i)=(x(i)-1)^2;
    if i>1
        S2(i-1)=x(i)*x(i-1);
    end
end
f=sum(S1)-sum(S2);
end