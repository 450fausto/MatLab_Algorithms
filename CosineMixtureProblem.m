function y=CosineMixtureProblem(x)
% Cosine Mixture Problem (CM) n=2 and n=4
% The number of local minima is one
% LimInf=[-1 -1]; LimSup=[1 1];
% Fojmin= -0.2 and -0.4 
% x*=[0 0] and [0 0 0 0];
S1=sum(cos(5*pi*x)); S2=sum(x.^2);
y=-0.1*S1+S2;
end