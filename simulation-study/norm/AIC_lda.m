function num=AIC_lda(U,mu0,mu1,lam,eigfun,arg,X,rang)
% select the number of cut-off level 
% for functional classification

ml=rang(1);mr=rang(2);
n=length(U);p_h=mean(U);
aic=zeros(1,mr-ml+1);
for m=ml:mr
    c1 = sum(lam(1:m).^(-1).*((mu0(1:m)).^2-(mu1(1:m)).^2))/2;
    c2 = (lam(1:m).^(-1).* (mu1(1:m)-mu0(1:m))) * eigfun(:,1:m)';
    eta = log(p_h/(1-p_h))+c1+trapz(arg, (X.*repmat(c2,n,1))');
    aic(m-ml+1)=-2*sum(U.*eta-log(1+exp(eta)))+2*(m+1);
end
[~,b]=min(aic); 
num=b+ml-1;