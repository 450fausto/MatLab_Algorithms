function f=Schaffer1(x)
% Schaffer 1 Problem (SF1) n=2
% LimInf=[-100 -100 ... -100]; LimSup=[100 100 ... 100];
% Fojmin=0
% x*=[0 0]
f=0.5+(((sin(norm(x)))^2)-0.5)/((1+0.001*(x(1)^2+x(2)^2 ))^2);
end