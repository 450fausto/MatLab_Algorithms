function f=Paviani(x)
% Paviani Problem (PP) n=10
% The number of local minima is unknown
% LimInf=[2 2 ... 2]; LimSup=[10 10 ... 10];
% Fojmin~=-45.778
% xi*~=9.351
S1=sum((log(x-2)).^2+(log(10-x)).^2);
P=prod(x);
f=S1-P^(1/5);
end