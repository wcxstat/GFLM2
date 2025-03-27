function xc = bisect(Y,U,rX,piX,family,a,b1,tol)

fb=1000*eqfun(Y,U,rX,piX,family,b1);
k=1;
if fb>0
    while fb>0
        k=k+1;
        fb=1000*eqfun(Y,U,rX,piX,family,k*b1);
    end
end
b=k*b1;

ind = b-a;
while ind > tol
    xx = (a+b)/2;
    fa=1000*eqfun(Y,U,rX,piX,family,a);
    fxx=1000*eqfun(Y,U,rX,piX,family,xx);
    if fa*fxx < 0
        b = xx;
    else
        a = xx;
    end
    ind = b - a;
end
xc = xx;