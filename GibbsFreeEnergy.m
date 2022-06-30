function y=GibbsFreeEnergy(k)
% White, W. B., Johnson, S. M., & Dantzig, G. B. (1958). Chemical 
% Equilibrium in Complex Mixtures. The Journal of Chemical Physics, 28(5), 
% 751–755. doi:10.1063/1.1744264 
% k=[0.0406681, 0.1477303, 0.7831534, 0.0014142, 0.4852466, 0.0006932, 0.0179473]
% LimInf=0
% LimSup=1
% Fobj=-47.761090859
n=k;
P=exp(3.932);
n(10)=2-n(1)-2*(n(2)+n(3))-n(6);
n(8)=1-n(4)-2*n(5)-n(6);
n(9)=0.5*(1-n(10)-n(3)-n(7)-n(8));
c=[-10.021, -21.096, -37.986, -9.846, -28.653, -18.918, -14.64, -28.032, -30.594, -26.111];
I=sum((log(n/sum(n))+log(P)+c).*n);
if imag(I)~=0
    y=abs(imag(I))+abs(real(I));
else
    y=I;
end
end