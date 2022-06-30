function [Solution,Convergence]=DEPrice(CostFunction,LimInf,LimSup,NumPop,MaxIter)
%Storn R. and Price K. (1997). Differential evolution: A simple and 
%efficient heuristic for global optimization over continuous spaces. 
%Journal of Global Optimization 11, 341–359.
%% At first...
if size(LimInf)~=size(LimSup)
    error('LimInf and LimSup must be the same size')
end
if size(LimSup,1)>1
    error('Both LimInf and LimSup must be a row vector')
end
if MaxIter<=0
    error('MaxIter must be greater than zero')
end
%% Control parameters
NumVar=length(LimSup);
F=0.5;
Cr=0.2;
%% Empty matrices
Fit=NaN(NumPop,1);
Convergence=zeros(MaxIter,1);
%% Poblation
X=repmat(LimInf,NumPop,1)+(repmat(LimSup,NumPop,1)-repmat(LimInf,NumPop,1)).*rand(NumPop,NumVar);
for i=1:NumPop
    Fit(i)=CostFunction(X(i,:));
end
[Fbest,b]=min(Fit);
Convergence(1)=Fbest;
iter=1;
Xbest=X(b,:);
%% Iterations
while iter<MaxIter
    for i=1:NumPop
        A=randperm(NumPop);
        A(A==i)=[];
        r1=A(1);
        r2=A(2);
        r3=A(3);
        RndInt=(randi(NumVar,[1,NumVar]))==(1:NumVar);
        RndCr=rand(1,NumVar)<Cr;
        n1=find(max(RndInt,RndCr)==1); %F=rand;
        y=X(r1,:)+F*(X(r2,:)-X(r3,:));
        z=X(i,:);
        z(n1)=y(n1);
        if prod(z<LimSup)==0 || prod(z>LimInf)==0
            z=LimInf+(LimSup-LimInf).*rand(1,NumVar);
        end
        FunZ=CostFunction(z);
        if FunZ<=Fit(i)
            X(i,:)=z;
            Fit(i)=FunZ;
            if FunZ<Fbest
                Xbest=z;
                Fbest=FunZ;
            end
        end
    end
    iter=iter+1;
    Convergence(iter)=Fbest;
end
Solution=[Xbest,Fbest];
end
