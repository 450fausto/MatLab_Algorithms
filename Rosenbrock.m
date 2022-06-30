function f=Rosenbrock(x)
% Rosenbrock Problem (RB) n=10
% LimInf=[-30 -30 ... -30]; LimSup=[30 30 ... 30];
% Fojmin=0
% x*=[1 1 ... 1]
n=length(x);
for i=1:(n-1)
    S(i)=100*(x(i+1)-x(i)^2)^2+(x(i)-1)^2;
end
f=sum(S);
end