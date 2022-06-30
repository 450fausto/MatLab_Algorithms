function y=Branin(x)
% Branin Problem (BP) n=2
% The number of local minima is three
% LimInf=[-5 0]; LimSup=[10 15];
% Fojmin= 5/(4*pi)
% x*=[-pi 12.275…],[pi 2.275…],[3*pi 2.475…];
a=1; b=5.1/(4*pi^2 ); c=5/pi; d=6; g=10; h=1/(8*pi);
y=a*(x(2)-b*x(1)^2+c*x(1)-d)^2+g*(1-h)*cos(x(1))+g;
end
