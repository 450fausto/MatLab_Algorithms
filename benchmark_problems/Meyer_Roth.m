function f=Meyer_Roth(x)
% Meyer and Roth Problem (MR) n=3
% The number of local minima is about one
% LimInf=[-10 -10 -10]; LimSup=[10 10 10];
% Fojmin=0.4 X 10^(-4)
% x*=[3.13 15.16 0.78]    <----- La solución está fuera de los límites
t=[1 2  1  2 0.1]; v=[1 1  2  2  0]; y=[0.126 0.219 0.076 0.126 0.186];
for i=1:5
f(i)=((x(1)*x(3)*t(i))/((1+x(1)*t(i)+x(2)*v(i)))-y(i))^2;
end
f=sum(f);
end
