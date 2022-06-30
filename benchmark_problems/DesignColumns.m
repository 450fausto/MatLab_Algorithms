function y=DesignColumns(x)
% Dinkoff, B., Levine, M., & Luus, R. (1979). Optimum Linear Tapering in 
% the Design of Columns. Journal of Applied Mechanics, 46(4), 956. 
% doi:10.1115/1.3424690 
% P, kN  L, m   E, GPa   omega     xi
% 10      2.5    10        1       1/12
% 
% a          Wo, m     W1, m     V, m3
% 0.02092    0.03606   0.06221   0.006178
P=10000; L=2.5; E=1e10; omega=1;
xi=1/12;
a=x(1); 
%W1=1.1849*sqrt(sqrt(P*L*L/(xi*pi*pi*E)));
W1=x(2);
k=sqrt(P/(xi*E));
h=k/(a*W1)-k/(a*(W1-a*L/2))+pi-atan(k/(a*W1));
V=omega*L*(W1^2-L*W1*a/2+((L*a)^2)/12)
y=V+1000*abs(h);
end
