function f=Shekel5(x)
% Shekel 5 Problem (S5) n=4
% LimInf=[0 0 0 0]; LimSup=[10 10 10 10];
% Fojmin~=-10.1499
% x*=[4 4 4 4]
n=length(x);
a=[4 4 4 4; 1 1 1 1; 8 8 8 8; 6 6 6 6; 3 7 3 7; 2 9 2 9; 5 5 3 3;...
    8 1 8 1; 6 2 6 2; 7 3.6 7 3.6];
c=[0.1 0.2 0.2 0.4 0.4 0.6 0.3 0.7 0.5 0.5]; S=zeros(1,n); S1=zeros(1,5);
for i=1:5
    for j=1:n
        S(j)=(x(j)-a(i,j))^2;
    end
    S1(i)=1/(sum(S)+c(i));
end
f=-sum(S1);
end
