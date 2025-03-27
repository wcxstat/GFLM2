function [res, it]=newtoneq(Y,U,rX,piX,family,x,precis)
% family=0: Gaussian with variance 1
% family=1: Logistic regression
% family=2: Poisson regression

it=0; x0=x;
if family==0
    dpsi1=rX;
    dpsi2=x0*rX;
    ddpsi1=ones(1,length(U));
    ddpsi2=ones(1,length(U));
    dddpsi=zeros(1,length(U));
elseif family==1
    dpsi1=expfun(rX);
    dpsi2=expfun(x0*rX);
    ddpsi1=expfun(rX).*(1-expfun(rX));
    ddpsi2=expfun(x0*rX).*(1-expfun(x0*rX));
    dddpsi=expfun(x0*rX).*(1-expfun(x0*rX)).*(1-2*expfun(x0*rX));
elseif family==2
    dpsi1=exp(rX);
    dpsi2=exp(x0*rX);
    ddpsi1=exp(rX);
    ddpsi2=exp(x0*rX);
    dddpsi=exp(x0*rX);
end
S1=ddpsi2.*piX;
S2=ddpsi1.*(1-piX);
eta=(1-U).*dpsi1+U.*dpsi2;
func=sum((S2.*U-(1-U).*S1*x0)./(S2+S1*x0^2).*(rX.*(Y-eta)));
dfunc1=sum(((1-U).*(S1*x0./(S2+S1*x0^2)).^2-...
    (S1.*S2)./((S2+S1*x0^2).^2).*(1-U+2*U*x0)-...
    dddpsi.* piX.*rX.*S2./((S2+S1*x0^2).^2).*(1-U+U*x0)*x0).*(rX.*(Y-eta)));
dfunc2=sum(S2./(S2+S1*x0^2).*U.*(rX.^2).*ddpsi2);
dfunc=dfunc1-dfunc2;
d=func/dfunc;
while abs(d)>precis
    x1=x0-d;
    it=it+1;
    x0=x1;
    if family==0
        dpsi1=rX;
        dpsi2=x0*rX;
        ddpsi1=ones(1,length(U));
        ddpsi2=ones(1,length(U));
        dddpsi=zeros(1,length(U));
    elseif family==1
        dpsi1=expfun(rX);
        dpsi2=expfun(x0*rX);
        ddpsi1=expfun(rX).*(1-expfun(rX));
        ddpsi2=expfun(x0*rX).*(1-expfun(x0*rX));
        dddpsi=expfun(x0*rX).*(1-expfun(x0*rX)).*(1-2*expfun(x0*rX));
    elseif family==2
        dpsi1=exp(rX);
        dpsi2=exp(x0*rX);
        ddpsi1=exp(rX);
        ddpsi2=exp(x0*rX);
        dddpsi=exp(x0*rX);
    end
    S1=ddpsi2.*piX;
    S2=ddpsi1.*(1-piX);
    eta=(1-U).*dpsi1+U.*dpsi2;
    func=sum((S2.*U-(1-U).*S1*x0)./(S2+S1*x0^2).*(rX.*(Y-eta)));
    dfunc1=sum(((1-U).*(S1*x0./(S2+S1*x0^2)).^2-...
        (S1.*S2)./((S2+S1*x0^2).^2).*(1-U+2*U*x0)-...
        dddpsi.* piX.*rX.*S2./((S2+S1*x0^2).^2).*(1-U+U*x0)*x0).*(rX.*(Y-eta)));
    dfunc2=sum(S2./(S2+S1*x0^2).*U.*(rX.^2).*ddpsi2);
    dfunc=dfunc1-dfunc2;
    d=func/dfunc;
end
res=x0;