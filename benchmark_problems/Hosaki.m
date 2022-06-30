function f=Hosaki(x)
% Hosaki Problem (HSK) n=2
% The number of local minima is two
% LimInf=[0 0]; LimSup=[5 6];
% Fojmin~=-2.3458 
% x*=[4 2]
f=(1-8*x(1)+7*x(1)^2-(7/3)*x(1)^3+(1/4)*x(1)^4)*(x(2)^2)*exp(-x(2));
end
