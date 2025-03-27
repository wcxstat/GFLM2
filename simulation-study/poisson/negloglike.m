function [value]=negloglike(Y,U,xi0,xi1,theta,m,param)

b00=param(1);
b10=param(2);
g=param(3:(m+2))';
com0=b00+xi0(:,1:m)*g;
com1=theta*(b10+xi1(:,1:m)*g);

value=-sum(Y(U==0).*com0'-(exp(com0')-1))-...
    sum(Y(U==1).*com1'-(exp(com1')-1));