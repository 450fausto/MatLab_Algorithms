function f=Sinusoidal(x)
% Sinusoidal Problem (SIN) n=10, 20
% LimInf=[0 0 ... 0]; LimSup=[180 180 ... 180];
% Fojmin~=-(A+1)
% x*=[90+z 90+z ... 90+z]
A=2.5; B=5; z=30; z=(pi/180)*z;
x=(pi/180)*x;
% n=length(x);       Number of local minima
% for i=1:(n/2)
%     Sol(i)=(factorial(n)/(factorial(n-2*i)*factorial(2*i)))*(3^(n-2*i))*(2^(2*i));
% end
% Sol=sum(Sol)
f=-(A*prod(sin(x-z))+prod(sin(B*(x-z))));
end
