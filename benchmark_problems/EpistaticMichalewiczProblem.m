function f=EpistaticMichalewiczProblem(x)
% Epistatic Michalewicz Problem (EM) n=5 and n=10
% The number of local minima is not known
% LimInf=[0 0 ... 0]; LimSup=[pi pi ... pi];
% Fojmin=-4.687658 and -9.660152
% x*=(2.693, 0.259, 2.074, 1.023, 1.720) and 
%  (2.693, 0.259, 2.074, 1.023, 2.275, 0.500, 2.138, 0.794, 2.219, 0.533)
Q=pi/6; m=10; n=length(x);
for i=1:n
    if rem(i,2)==1 && i<n
        y(i)=x(i)*cos(Q)-x(i+1)*sin(Q);
    elseif rem(i,2)==0 
        y(i)=x(i-1)*sin(Q)+x(i)*cos(Q);
     elseif i==n && rem(i,2)==1
         y(i)=x(i);
    end
    f(i)=sin(y(i))*(sin((i*y(i)^2)/pi))^(2*m);
end
f=-sum(f);
end

