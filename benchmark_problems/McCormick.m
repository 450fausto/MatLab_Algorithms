function f=McCormick(x)
% McCormick Problem (MC) n=2
% The number of local minima is about one
% LimInf=[-1.5 -3]; LimSup=[4 3];
% Fojmin=-1.9133
% x*=[-0.547 -1.547]
f=sin(x(1)+x(2))+(x(1)-x(2))^2-(3/2)*x(1)+(5/2)*x(2)+1;
end
