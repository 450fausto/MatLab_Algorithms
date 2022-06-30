function f=Hartman3(x)
% Hartman 3 Problem (H3) n=3
% The number of local minima is four
% LimInf=[0 0 0]; LimSup=[1 1 1];
% Fojmin?-3.862782
% x*?[0.114614,0.555649,0.852547]
a=[3 10 30; 0.1 10 35; 3 10 30; 0.1 10 35];
p=[0.3689 0.117 0.2673; 0.4699 0.4387 0.747; 0.1091 0.8732 0.5547; ... 
    0.03815 0.5743 0.8828];
c=[1 1.2 3 3.2]; 
for i=1:4
    for j=1:3
        f1(j)=a(i,j)*(x(j)-p(i,j))^2;
    end
    s1=-sum(f1);
    f2(i)=c(i)*exp(s1);
end
f=-sum(f2);
end
