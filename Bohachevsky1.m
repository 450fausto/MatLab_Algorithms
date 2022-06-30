function y=Bohachevsky1(x)
% Bohachevsky 1 Problem (BF1) n=2
% The number of local minima is unknown
% LimInf=[-50 -50]; LimSup=[50 50];
% Fojmin= 0
% x*=[0 0];
y=x(1)^2+2*x(2)^2-0.3*cos(3*pi*x(1))-0.4*cos(4*pi*x(2))+0.7;
end