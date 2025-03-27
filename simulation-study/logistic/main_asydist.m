p=0.6;
%n=2000;
nvec=[200,600,2000];
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


mu0_c=zeros(1,M);mu1_c=1*(-1).^(2:(M+1)).*(1:M).^(-3);
mu=[mu0_c;mu1_c]*Phi;

asynorm_mat=zeros(length(nvec),1000);
for nn=1:length(nvec)
n=nvec(nn);
theta_mat=zeros(3,1000);
mise_mat=zeros(4,1000);
asynorm=zeros(1,1000);
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

Y=zeros(1,n);
prob=1./(1+exp(-(1-U+U*theta).*Int));
for i=1:n
    Y(i)=binornd(1,prob(i));
end

%%% Common FPCA %%%
p_h=mean(U);
X0=X(U==0,:); X1=X(U==1,:);
mu0=mean(X0); mu1=mean(X1);
X_cen0=X0-repmat(mu0,n0,1); % center
X_cen1=X1-repmat(mu1,n1,1);
Xcov0=(X_cen0)'*X_cen0/n0; % covaraince matrix
Xcov1=(X_cen1)'*X_cen1/n1;
Xcov_p=(1-p_h)*Xcov0+p_h*Xcov1; % joint
[lam,eigf]=inteq(Xcov_p,arg,15); % eigvalue and eigfunction
xi0=zeros(n0,15);xi1=zeros(n1,15); % FPC score
for j=1:15
    xi0(:,j)=trapz(arg,(X_cen0.*repmat(eigf(:,j)',n0,1))');
    xi1(:,j)=trapz(arg,(X_cen1.*repmat(eigf(:,j)',n1,1))');
end

%%% Estimation of r(X_i) %%%
% m1 = 3;
Y0=Y(U==0); Y1=Y(U==1);
m1=AIC(Y0,xi0,1,[1,14]);
%m1=2;
bh_coef=glmfit(xi0(:,1:m1),Y0','binomial');
b_h=bh_coef(2:end)'*eigf(:,1:m1)';
inter_h=bh_coef(1)-trapz(arg,mu0.*b_h);
Int_h=inter_h+trapz(arg,(X.*repmat(b_h,n,1))');

%%%  Estimation of pi(X_i) by LDA  %%%
mu0_coef=trapz(arg, eigf.*repmat(mu0',1,15));
mu1_coef=trapz(arg, eigf.*repmat(mu1',1,15));
m2=AIC_lda(U,mu0_coef,mu1_coef,lam,eigf,arg,X,[1,14]);
% m2 = 4;
c1=sum(lam(1:m2).^(-1).*((mu0_coef(1:m2)).^2-(mu1_coef(1:m2)).^2))/2;
c2=(lam(1:m2).^(-1).* (mu1_coef(1:m2)-mu0_coef(1:m2))) * eigf(:,1:m2)';
Q_lda=log(p_h/(1-p_h))+c1+trapz(arg, (X.*repmat(c2,n,1))');
pi_lda=expfun(Q_lda); % a 1-by-n vector

%%%  Estimation of pi(X_i) using true value  %%%
c1=sum(lamsq(1,:).^(-2).*(mu0_c.^2-mu1_c.^2)/2);
c2=lamsq(1,:).^(-2).*(mu1_c-mu0_c) * Phi;
Q_true=log(p/(1-p))+c1+trapz(arg, (X.*repmat(c2,n,1))');
pi_true=expfun(Q_true); % a 1-by-n vector

theta_h1=bisect(Y,U,Int_h,pi_true,1,0,5,0.001); % true
theta_h2=bisect(Y,U,Int_h,pi_lda,1,0,5,0.001); % lda
theta_h3=bisect(Y,U,Int_h,p_h,1,0,5,0.001); % ind

bh_coef1=glmfit(xi1(:,1:m1),Y1','binomial');
% inital value
x_in1=[bh_coef(1),bh_coef1(1)/theta_h1,...
    (1-p_h)*bh_coef(2:end)'+p_h*bh_coef1(2:end)'/theta_h1];
x_in2=[bh_coef(1),bh_coef1(1)/theta_h2,...
    (1-p_h)*bh_coef(2:end)'+p_h*bh_coef1(2:end)'/theta_h2];
x_in3=[bh_coef(1),bh_coef1(1)/theta_h3,...
    (1-p_h)*bh_coef(2:end)'+p_h*bh_coef1(2:end)'/theta_h3];

options=optimset('MaxFunEvals', 1000, 'MaxIter', 1000,'TolX',1e-10,'TolFun',1e-10, 'Display', 'off' );
bh_coef1=fmincon(@(param)negloglike(Y,U,xi0,xi1,theta_h1,m1,param), x_in1, [],[], [],[], [], [],[], options);
bh_coef2=fmincon(@(param)negloglike(Y,U,xi0,xi1,theta_h2,m1,param), x_in2, [],[], [],[], [], [],[], options);
bh_coef3=fmincon(@(param)negloglike(Y,U,xi0,xi1,theta_h3,m1,param), x_in3, [],[], [],[], [], [],[], options);
br_h1=bh_coef1(3:end)*eigf(:,1:m1)';
br_h2=bh_coef2(3:end)*eigf(:,1:m1)';
br_h3=bh_coef3(3:end)*eigf(:,1:m1)';

%Intr_h=bh_coef2(1)+trapz(arg,((X-mu0).*repmat(br_h2,n,1))');
Ibd=ibound(Int_h,pi_lda,1,theta_h2);
asynorm(num)=sqrt(n*Ibd)*(theta_h2-theta);

theta_mat(:,num)=[theta_h1;theta_h2;theta_h3];
mise_mat(:,num)=[trapz(arg,(b_h-b).^2);trapz(arg,(br_h1-b).^2);...
    trapz(arg,(br_h2-b).^2);trapz(arg,(br_h3-b).^2)];
end
asynorm_mat(nn,:)=asynorm;
end

ksdensity(asynorm_mat(1,:))
hold on
ksdensity(asynorm_mat(2,:))
hold on
ksdensity(asynorm_mat(3,:))
y=normpdf(sort(asynorm_mat(3,:)),0,1);
plot(sort(asynorm_mat(3,:)),y,'--r')

mean(abs(asynorm_mat)<=norminv(0.975),2)

save('asynorm_logit_alpha2.0.txt','asynorm_mat','-ascii')