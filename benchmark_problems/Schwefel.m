function f=Schwefel(x)
% Schwefel Problem (SWF)  n=10
% LimInf=[-500 -500 ... -500]; LimSup=[500 500 ... 500];
% Fojmin~=-418.9829n
% x*=[420.97 420.97 ... 420.97]
f=-sum(x.*sin(sqrt(abs(x))));
end
