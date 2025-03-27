p=0.6;
n=600;
theta=1.5;
alpha=2.0;
arg=0:0.01:1;
L=length(arg);

M=50;
Phi=ones(M,L);
for i=2:M
    Phi(i,:)=sqrt(2)*cos((i-1)*pi*arg);
end
lamsq=[((-1).^(2:(M+1))).*((1:M).^(-alpha/2));
    ((-1).^(2:(M+1))).*((1:M).^(-alpha/2))];
b_coef=1*(-1).^(2:(M+1)).*((1:M).^(-3));
b=b_coef*Phi; % true slope function

gamma_vec=[100,3,1.5];
res_b=zeros(1,6);
for rep=1:3
fprintf(' Case %d\n',rep);
if rep~=1
mu0_c=zeros(1,M);mu1_c=1*(-1).^(2:(M+1)).*(1:M).^(-gamma_vec(rep));
else
mu0_c=zeros(1,M);mu1_c=zeros(1,M);
end
mu=[mu0_c;mu1_c]*Phi;

mise_mat=zeros(1,1000);
for num=1:1000
X=zeros(n,L);
U=binornd(1,p,1,n);
n0=sum(1-U);n1=sum(U);
Int=zeros(1,n); 
for i=1:n
    id=U(i)+1;
    X(i,:)=mu(id,:)+(normrnd(0,1,1,M).*lamsq(id,:))*Phi;
    Int(i)=trapz(arg,X(i,:).*b);
end

Y=(1-U+U*theta).*Int+normrnd(0,1,1,n);

mu_hat=mean(X);
X_cen=X-repmat(mu_hat,n,1); % center
Xcov=(X_cen)'*X_cen/n; % covaraince matrix
[lam,eigf]=inteq(Xcov,arg,15); % eigvalue and eigfunction
xi=zeros(n,15);% FPC score
for j=1:15
    xi(:,j)=trapz(arg,(X_cen.*repmat(eigf(:,j)',n,1))');
end
m1=AIC(Y,xi,0,[1,14]);
%m1=2;
bh_coef=glmfit(xi(:,1:m1),Y,'normal');
b_h=bh_coef(2:end)'*eigf(:,1:m1)';

mise_mat(num)=trapz(arg,(b_h-b).^2);
end
res_b(((rep-1)*2+1):(rep*2))=[mean(mise_mat);std(mise_mat)];
end
res_b