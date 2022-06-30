function f=Kowalik(x)
% Kowalik Problem (KL) n=4
% The number of local minima is one
% LimInf=[0 0 0 0]; LimSup=[0.42 0.42 0.42 0.42];
% Fojmin~=3.0748 X 10^(-4)
% x*~=[0.192 0.190 0.123 0.135]
a=[0.1957 0.1947 0.1735 0.16 0.0844 0.0627 0.0456 0.0342 0.0323 0.0235 0.0246];
b=[0.25 0.50 1.0 2.0 4.0 6.0 8.0 10.0 12.0 14.0 16.0];
for i=1:11
    f(i)=(a(i)-(x(1)*(1+x(2)*b(i)))/(1+x(3)*b(i)+x(4)*b(i)^2))^2;
end
f=sum(f);
end