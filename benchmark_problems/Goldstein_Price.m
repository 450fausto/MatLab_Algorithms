function f=Goldstein_Price(x)
% Goldstein and Price (GP) n=2
% The number of local minima is one
% LimInf=[-2 -2]; LimSup=[2 2];
% Fojmin=3
% x*=[0 -1]
f1=(1+((x(1)+x(2)+1)^2)*...
    (19-14*x(1)+3*x(1)^2-14*x(2)+6*x(1)*x(2)+3*x(2)^2));
f2=(30+((2*x(1)-3*x(2))^2)*...
    (18-32*x(1)+12*x(1)^2+48*x(2)-36*x(1)*x(2)+27*x(2)^2));
f=f1*f2;
end
