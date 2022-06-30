function f=HelicalValley(x)
% Helical Valley Problem (HV) n=3
% The number of local minima is one
% LimInf=[-10 -10 -10]; LimSup=[10 10 10];
% Fojmin=0
% x*=[1 0 0]
if x(1)>=0
    Q=(1/(2*pi))*atan(x(2)/x(1));
elseif x(1)<0
    Q=(1/(2*pi))*atan(x(2)/x(1))+1/2;
end
f=100*((x(2)-10*Q)^2+(sqrt(x(2)^2+x(1)^2)-1)^2)+x(3)^2;
end