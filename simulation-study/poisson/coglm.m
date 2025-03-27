function [coef,it]=coglm(Y,U,xi0,xi1,theta,family,m,x,precis)
% family=0: Gaussian with variance 1
% family=1: Logistic regression
% family=2: Poisson regression
% Input:
%    xi0 correponds to U=0
%    xi1 correponds to U=1
%    x is the initial value
%    precis is the required accuracy

% x=zeros(1,m1+2);m=m1;
% family=2;theta=1.5;
it=0; x0=x;
[res,hess]=ddev(Y,U,xi0,xi1,theta,family,m,x0);
d=res/hess;
while norm(d)>precis
    x1=x0-d;
    it=it+1;
    x0=x1;
    [res,hess]=ddev(Y,U,xi0,xi1,theta,family,m,x0);
    d=res/hess;
end
coef=x0;
