function f=Rastrigin(x)
% Rastrigin Problem (RG) n=10
% LimInf=[-5.12 -5.12 ... -5.12]; LimSup=[5.12 5.12 ... 5.12];
% Fojmin=0
% x*=[0 0 ... 0]
n=length(x);
f=10*n+sum(x.^2-10*cos(2*pi*x));
end
