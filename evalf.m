function y=evalf(X,Y,C,t,w,z)
y=0.;
m=size(X,1);
n=size(X,2);

for i=1:m
    if(i <= n)
        y = y+t/2*w(i)*w(i);
    end
    y = y+ t*C*z(i);
    %y=y+t*1./2*w(i)^2 + t*C*z(i);
    vv = Y(i)*X(i,:)*w + z(i) - 1;
    if(z(i) <= 0 || vv <= 0)
        y=inf;
        break;
    end
    y=y-log(z(i));
    y=y-log(vv);
end
