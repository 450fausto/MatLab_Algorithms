function f=Griewank(x)
% Griewank Problem (GW) n=10
% The number of local minima is unknown
% LimInf=[-600 -600 ... -600]; LimSup=[600 600... 600];
% Fojmin=0
% x*=[0 0 ... 0]
i=1:length(x);
f=1+(1/4000)*sum(x.^2)-prod(cos(x./sqrt(i)));
end