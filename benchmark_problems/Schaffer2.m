function f=Schaffer2(x)
% Schaffer 2 Problem (SF1) n=2
% LimInf=[-100 -100 ... -100]; LimSup=[100 100 ... 100];
% Fojmin=0
% x*=[0 0]
f=((x(1)^2+x(2)^2)^0.25)*(sin(50*(x(1)^2+x(2)^2)^0.1)^2+1);
end
