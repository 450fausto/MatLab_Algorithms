function f=PowellQuadratic(x)
% Powell’s Quadratic Problem (PWQ) n=4
% LimInf=[-10 -10 -10 -10]; LimSup=[10 10 10 10];
% Fojmin=0
% x*=[0 0 0 0]
f=(x(1)+10*x(1))^2+5*(x(3)-x(4))^2+(x(2)-2*x(3))^4+10*(x(1)-x(4))^4;
end
