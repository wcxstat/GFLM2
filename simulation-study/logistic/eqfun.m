function res=eqfun(Y,U,rX,piX,family,rang)
% family=0: Gaussian with variance 1
% family=1: Logistic regression
% family=2: Poisson regression

L=length(rang);
n=length(Y);
result=zeros(1,L);
for i=1:L
    x0=rang(i);
    if family==0
        dpsi1=rX;
        dpsi2=x0*rX;
        ddpsi1=ones(1,n);
        ddpsi2=ones(1,n);
    elseif family==1
        dpsi1=expfun(rX);
        dpsi2=expfun(x0*rX);
        ddpsi1=expfun(rX).*(1-expfun(rX));
        ddpsi2=expfun(x0*rX).*(1-expfun(x0*rX));
    elseif family==2
        dpsi1=exp(rX);
        dpsi2=exp(x0*rX);
        ddpsi1=exp(rX);
        ddpsi2=exp(x0*rX);
    end
    S1=ddpsi2.*piX;
    S2=ddpsi1.*(1-piX);
    ep=Y-(1-U).*dpsi1-U.*dpsi2;
    result(i)=sum((S2.*U-(1-U).*S1*x0)./(S2+S1*x0^2).*rX.*ep);
end
res=result;