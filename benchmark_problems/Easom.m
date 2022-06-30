function y=Easom(x)
% Easom Problem (EP) n=2
% The number of local minima is one
% LimInf=[-10 -10]; LimSup=[10 10];
% Fojmin=-1
% x*=[pi pi] 
y=-cos(x(1))*cos(x(2))*exp(-(x(1)-pi)^2-(x(2)-pi)^2);
end
