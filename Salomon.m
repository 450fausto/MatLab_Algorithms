function f=Salomon(x)
% Salomon Problem (SAL) n=5, 10
% LimInf=[-100 -100 ... -100]; LimSup=[100 100 ... 100];
% Fojmin=0
% x*=[0 0 ... 0]
f=1-cos(2*pi*norm(x))+norm(x)/10;
end