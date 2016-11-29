function [w,z] = central(X,Y,C,t,w0,z0)

m=size(w0,1);
epsilon = 0.0001;
w=w0;
z=z0;
hessi=zeros(2*m,2*m);
grad=zeros(2*m,1);

while(true)
    precalc=zeros(m,1);
    
    for i = 1:m
        precalc(i) = Y(i)*(X(i,:)*w) + z(i) - 1;
        %precalc(i) = precalc(i) * precalc(i);
        %X(i,:) = X(i,:)/precalc(i);
    end
    
	for i = 1:(m)
        grad(i)=t*w(i);
        grad(i+m)=-1/z(i) +t*C - 1/precalc(i);
		for j=1:(m)
			%premier cas : wi wj
			%hessi(i,j) = X(:,i)'*X(:,j) + (i==j)*t;
            hessi(i,j)=(i==j)*t;
            for k=1:m
                hessi(i,j) = hessi(i,j) + X(k,i)*X(k,j)/(precalc(k)^2);
            end
            %deuxi?me cas : zi zj
            hessi(i+m, j+m) = (i==j)*(1/(z(i)*z(i)) + 1. /(precalc(i) ^ 2));
            %for p=1:m
            %    hessi(i+m, j+m) = hessi(i+m, j+m) + 1./(precalc(p)*precalc(p));
            %end
            
            %troisi?me cas : wi zj
            hessi(i,j+m) = Y(j)*X(j,i)/(precalc(j)^2);
            hessi(i+m, j) = hessi(i, j+m);
            
            %calc du grad i
            grad(i) = grad(i) - Y(j)*X(j,i)/precalc(j);
            %grad(i+m) = grad(i+m) - 1/precalc(j);
		end
    end
    
    for i = 1:m
        %X(i,:) = X(i,:)*precalc(i);
    end
    
    invhess = inv(hessi);
    %x  = [w; z];
    dx = -invhess * grad;
    lambda = -grad'*dx;
    if(lambda <= epsilon*2)
        break;
    end
    sz=1;
    alpha=1/4;
    beta=1/2;
    while (evalf(X,Y,C,t,w+sz*dx(1:m),z+sz*dx((m+1):(2*m))) >= evalf(X,Y,C,t,w,z) + alpha*sz*grad'*dx)
        sz=beta*sz;
    end
    %x = x+sz*dx;
    w  = w+sz*dx(1:m);
    z  = z+sz*dx((m+1):(2*m));
end

end
