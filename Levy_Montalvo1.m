function f=Levy_Montalvo1(x)
% Levy and Montalvo 1 Problem (LM1) n=3
% The number of local minima is about 5^n
% LimInf=[-10 -10 -10]; LimSup=[10 10 10];
% Fojmin=0
% x*=[-1 -1 -1]
n=length(x);
for i=1:n
    y(i)=1+(1/4)*(x(i)+1);
end
for i=1:(n-1)
    S(i)=((y(i)-1)^2)*(1+10*(sin(pi*y(i+1)))^2);
end
S=sum(S);
f=(pi/n)*(10*(sin(pi*y(1)))^2)+S+(y(n)-1)^2;
end