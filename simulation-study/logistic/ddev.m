function [res,hess]=ddev(Y,U,xi0,xi1,theta,family,m,x0)
% xi0 and xi1 are n0-by-M and n1-by-M matrices

% Input:
%   x0: a 1-by-(2+m) vector
% Output:
%   res: a 1-by-(2+m) vector
%   hess: a (2+m)-by-(2+m) Hessian matrix

b00=x0(1);
b10=x0(2);
g=x0(3:(m+2))';
com0=b00+xi0(:,1:m)*g;
com1=theta*(b10+xi1(:,1:m)*g);
[dpsi0,ddpsi0]=linkfun(family,com0');
[dpsi1,ddpsi1]=linkfun(family,com1');

res=zeros(1,2+m);
ep0=Y(U==0)-dpsi0;
ep1=Y(U==1)-dpsi1;
res(1)=sum(ep0);
res(2)=sum(ep1)*theta;
res(3:end)=xi0(:,1:m)'*ep0'+theta*(xi1(:,1:m)'*ep1');

hess=zeros(2+m,2+m);
hess(1,1)=-sum(ddpsi0);
hess(1,3:end)=-xi0(:,1:m)'*ddpsi0';
hess(2,2)=-theta^2*sum(ddpsi1);
hess(2,3:end)=-theta^2*(xi1(:,1:m)'*ddpsi1');
hess(3:end,1)=hess(1,3:end)';
hess(3:end,2)=hess(2,3:end)';
for i=1:sum(1-U)
    hess(3:end,3:end)=hess(3:end,3:end)-ddpsi0(i)*xi0(i,1:m)'*xi0(i,1:m);
end
for i=1:sum(U)
    hess(3:end,3:end)=hess(3:end,3:end)-theta^2*ddpsi1(i)*xi1(i,1:m)'*xi1(i,1:m);
end
res=res/length(U);
hess=hess/length(U);