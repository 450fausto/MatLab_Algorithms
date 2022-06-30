function f=Exponential(x)
% Exponential Problem (EXP) n=10
% The number of local minima is one
% LimInf=[-1 -1 ... -1]; LimSup=[1 1 ... 1];
% Fojmin=-1
% x*=[0 0 ... 0]
f=-exp(-0.5*sum(x.^2));
end