function num=AIC(Y,X,family,rang)
% Y is a 1-by-n vector
% X is a n-by-p matrix
% range is 1-by-2 vector, i.e., range=[ml,mr]

n=length(Y);
ml=rang(1);mr=rang(2);
aic=zeros(1,mr-ml+1);
for m=ml:mr
    if family==0
        coef=glmfit(X(:,1:m),Y,'normal');
        eta=coef(1)+coef(2:end)'*X(:,1:m)';
        aic(m-ml+1)=-2*sum(Y.*eta-0.5*(eta.^2))+2*(m+1);
    elseif family==1
        coef=glmfit(X(:,1:m),Y','binomial');
        eta=coef(1)+coef(2:end)'*X(:,1:m)';
        aic(m-ml+1)=-2*sum(Y.*eta-log(1+exp(eta)))+log(n)*(m+1);
    else
        coef=glmfit(X(:,1:m),Y,'poisson');
        eta=coef(1)+coef(2:end)'*X(:,1:m)';
        aic(m-ml+1)=-2*sum(Y.*eta-exp(eta))+2*(m+1);
    end
end
[~,b]=min(aic); 
num=b+ml-1;