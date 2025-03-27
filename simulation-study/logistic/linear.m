n = 200;
alpha = 1.5;
arg = 0:0.01:1; L = length(arg);

M = 50;
Phi = ones(M,L);
for i = 2:M
    Phi(i,:) = sqrt(2)*cos((i-1)*pi*arg);
end

lamsq = ((-1).^(2:(M+1))) .* ((1:M).^(-alpha/2));
b_coef = 4*(-1).^(2:(M+1)) .* ((1:M).^(-3));
b = b_coef * Phi;

mise=zeros(1,1000);
for num = 1:1000
X = zeros(n,L);
xi_t=zeros(n,M);
Int = zeros(1,n); 
for i = 1:n
    xi_t(i,:)=normrnd(0,1,1,M) .* lamsq;
    X(i,:) = xi_t(i,:) * Phi;
    Int(i) = trapz(arg,X(i,:).*b);
end

Y=zeros(1,n);
prob=expfun(Int);
for i=1:n
    Y(i)=binornd(1,prob(i));
end

% mu = mean(X);
% X_cen = X - repmat(mu,n,1); % center
% Xcov = (X_cen)' * X_cen/n; % covaraince matrix
% [lam, eigf] = inteq(Xcov, arg, 15); % eigvalue and eigfunction
% xi = zeros(n, 15); % FPC score
% for j = 1:15
%     xi(:,j) = trapz(arg, (X_cen .* repmat(eigf(:,j)',n,1))');
% end
[lam, eigf, xi]=FPCA_bal(X, arg, 20);

% m1=AIC(Y,xi,1,[1,14]);
m1 = 3;
bh_coef=glmfit(xi(:,1:m1),Y','binomial');
b_h = bh_coef(2:end)'*eigf(:,1:m1)';
% b_h = [4,0.5]*eigf(:,1:m1)';
% inter_h = bh_coef(1)-trapz(arg,mu.*b_h);
% plot(arg,b_h)
% hold on
% plot(arg,b)
mise(num)=trapz(arg,(b_h-b).^2);
end
