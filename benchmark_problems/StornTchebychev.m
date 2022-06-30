function f=StornTchebychev(x)                                        
% Storn’s Tchebychev Problem (ST) n=9, 17                            |
% LimInf=[-128 -128 ... -128]; LimSup=[128 128 ... 128];             |
% Fojmin=0                                                          \|/ 
% x*=[128 0 -256 0 160 0 -32 0,1]  la dimensión de este vector es 16 V
% [32768 0 -1331072 0 21299 0 -180224 84480 0 -2154 0 2688 0 -128 0 1]
n=length(x); m=60*(n==9)+100*(n==17); k=0; 
d=72.661*(n==9)+10558.145*(n==17);
for j=1:(m+1)
    k=k+1;
    for i=1:n
        u(i)=((1.2)^(n-i))*x(i);
        v(i)=((-1.2)^(n-i))*x(i);
        s(i)=((2*(j-1)/m-1)^(n-i))*x(i);
    end
    u=sum(u); v=sum(v);
    p1=((u-d)^2)*(u<d); 
    p2=((v-d)^2)*(v<d); 
    w(j)=sum(s);
    p(j)=((w(j)-1)^2)*(w(j)>1)+((w(j)+1)^2)*(w(j)<-1);
end
p3=sum(p);
f=p1+p2+p3;
end
