function res=ibound(rX,piX,family,theta_t)
% family=0: Gaussian with variance 1
% family=1: Logistic regression
% family=2: Poisson regression

n=length(rX);
x0=theta_t;
if family==0
    ddpsi1=ones(1,n);
    ddpsi2=ones(1,n);
elseif family==1
    ddpsi1=exp(rX)./((1+exp(rX)).^2);
    ddpsi2=exp(x0*rX)./((1+exp(x0*rX)).^2);
elseif family==2
    ddpsi1=exp(rX);
    ddpsi2=exp(x0*rX);
end
S1=ddpsi2.*piX;
S2=ddpsi1.*(1-piX);
res=mean((S1.*S2.*(rX.^2))./(S2+S1*x0^2));