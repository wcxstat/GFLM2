function res=expfun(x)
% calculate res=exp(x)/(1+exp(x))

n=length(x);
res=zeros(1,n);
for i=1:n
    if x(i)>3
        res(i)=1/(1+exp(-x(i)));
    else
        res(i)=exp(x(i))/(1+exp(x(i)));
    end
end