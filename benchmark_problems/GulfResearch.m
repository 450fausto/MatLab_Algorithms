function f=GulfResearch(x)
% Gulf Research Problem (GRP) n=3
% The number of local minima is one
% LimInf=[0.1 0 0]; LimSup=[100 25.6 5];
% Fojmin=0
% x*=[50,25,1.5]
m=99;
for i=1:m
    u=25+(-50*log(i/100))^(1/1.5);
    f(i)=(exp(-((u-x(2))^x(3))/(x(1)))-i/100)^2;
end
f=sum(f);
end
