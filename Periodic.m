function f=Periodic(x)
% Periodic Problem (PRD) n=2
% The number of local minima is 49
% LimInf=[-10 -10]; LimSup=[10 10];
% Fojmin=0.9
% x*=[0 0]
f=1+(sin(x(1)))^2+(sin(x(2)))^2-(1/10)*exp(-x(1)^2-x(2)^2);
end