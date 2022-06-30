function f=Miele_Cantrell(x)
% Miele and Cantrell Problem (MCP) n=4
% The number of local minima is unknown
% LimInf=[-1 -1 -1 -1]; LimSup=[1 1 1 1];
% Fojmin=0
% x*=[0 1 1 1] 
f=(exp(x(1))-x(2))^4+100*(x(2)-x(3))^6+(tan(x(3)-x(4)))^4+x(1)^8;
end
