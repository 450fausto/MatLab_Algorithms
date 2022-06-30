function f=ModifiedLangerman(x)
% Modified Langerman Problem (ML) n=10
% The number of local minima is unknown
% LimInf=[0 0 ... 0]; LimSup=[10 10 ... 10];
% Fojmin=-0.965
% x*=[8.074 8.777 3.467 1.867 6.708 6.349 4.534 0.276 7.633 1.567]
n=length(x);
c=[0.806 0.517 0.1 0.908 0.965];
a=[9.681,0.667,4.783,9.095,3.517,9.325,6.544,0.211,5.122,2.02;9.4,2.041,...
    3.788,7.931,2.882,2.672,3.568,1.284,7.033,7.374;8.025,9.152,5.114,...
    7.621,4.564,4.711,2.996,6.126,0.734,4.982;2.196,0.415,5.649,6.979,...
    9.51,9.166,6.304,6.054,9.377,1.426;8.074,8.777,3.467,1.867,6.708,...
    6.349,4.534,0.276,7.633,1.567];
S1=zeros(1,n); d=zeros(1,5); S2=zeros(1,5);
for j=1:5
    for i=1:n
        S1(i)=(x(i)-a(j,i))^2;
    end
    d(j)=sum(S1);
    S2(j)=c(j)*cos(d(j)/pi)*exp(-pi*d(j));
end
f=-sum(S2);
end