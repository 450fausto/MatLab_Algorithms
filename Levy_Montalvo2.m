function f=Levy_Montalvo2(x)
% Levy and Montalvo 2 Problem (LM2) n=5, 10
% The number of local minima is about 15^n
% LimInf=[-5 -5 ... -5]; LimSup=[5 5 ... 5];
% Fojmin=0
% x*=[1 1 ... 1]
n=length(x);
for i=1:(n-1)
    S(i)=((x(i)-1)^2)*(1+(sin(3*pi*x(i+1)))^2);
end
S=sum(S);
f=(1/10)*((sin(3*pi*x(1)))^2+S+((x(n)-1)^2)*(1+(sin(2*pi*x(n)))^2));
end