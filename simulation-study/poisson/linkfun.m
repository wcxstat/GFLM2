function [dfun,ddfun]=linkfun(family,x0)
% family=0: Gaussian with variance 1
% family=1: Logistic regression
% family=2: Poisson regression

if family==0
    dfun=x0;
    ddfun=ones(1,length(x0));
elseif family==1
    dfun=expfun(x0);
    ddfun=expfun(x0).*(1-expfun(x0));
else
    dfun=exp(x0);
    ddfun=exp(x0);
end