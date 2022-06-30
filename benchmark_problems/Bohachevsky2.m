function y=Bohachevsky2(x)
% Bohachevsky 2 Problem (BF2) n=2
% The number of local minima is unknown
% LimInf=[-50 -50]; LimSup=[50 50];
% Fojmin= 0
% x*=[0 0];
y=x(1)^2+2*x(2)^2-0.3*cos(3*pi*x(1))*cos(4*pi*x(2))+0.3;
end
