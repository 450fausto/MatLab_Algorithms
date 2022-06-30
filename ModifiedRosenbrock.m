function f=ModifiedRosenbrock(x)
% Modified Rosenbrock Problem (MRP) n=2
% The number of local minima is two
% LimInf=[-5 -5]; LimSup=[5 5];
% Fojmin=0
% x*=[0.3412 0.1164] [1 1]
f=100*(x(2)-x(1)^2 )^2+(6.4*(x(2)-0.5)^2-x(1)-0.6)^2;
end