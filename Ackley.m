function y=Ackley(x)
% Ackley's problem (ACK) n=10
% The number of local minima is not known
% LimInf=[-30 -30 ... -30]; LimSup=[30 30 ... 30];
% Fojmin=0
% x*=[0 0 ... 0];
n=length(x);
S1=sum(x.^2);
S2=sum(cos(2*pi*x));
y=-20*exp(-0.02.*sqrt(S1/n))-exp(S2/n)+20+exp(1);
end