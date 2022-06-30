function f=MultiGaussian(x)
% Multi-Gaussian Problem (MGP) n=2
% The number of local minima is four
% LimInf=[-2 -2]; LimSup=[2 2];
% Fojmin~=-1.29695
% x*=[-0.01336 -0.01356]
a=[0.5 1.2 1 1 1.2]; b=[0 1 0 -0.5 0]; c=[0 0 -0.5 0 1]; 
d=[0.1 0.5 0.5 0.5 0.5];
S=zeros(1,5);
for i=1:5
    S(i)=a(i)*exp(-(((x(1)-b(i))^2+(x(2)-c(i))^2))/(d(i)^2));
end
f=-sum(S);
end
