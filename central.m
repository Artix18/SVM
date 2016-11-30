function [w,z] = central(X,Y,C,t,w0,z0)

n=size(X,2);
m=size(X,1);
epsilon = 0.0001;
w=w0;
z=z0;


NN = max(n,m);

while(true)
    precalc=zeros(m,1);
    hessi=zeros(n+m,n+m);
    grad=zeros(n+m,1);
    
    for i = 1:m
        precalc(i) = Y(i)*(X(i,:)*w) + z(i) - 1;
        %precalc(i) = precalc(i) * precalc(i);
        %X(i,:) = X(i,:)/precalc(i);
    end
    
	for i = 1:(n+m)
        if(i <= n)
            grad(i)=t*w(i);
        else
            grad(i)=-1/z(i-n) +t*C - 1/precalc(i-n);
        end
		for j=1:(n+m)
            %%% HESSIENNE %%%
			%premier cas : wi wj; i <= n et j <= n !!
			%hessi(i,j) = X(:,i)'*X(:,j) + (i==j)*t;
            if i <= n && j <= n
                hessi(i,j)=(i==j)*t;
                for k=1:m
                   hessi(i,j) = hessi(i,j) + X(k,i)*X(k,j)/(precalc(k)^2);
                end
            elseif i > n && j > n
            %deuxi?me cas : zi zj : i,j > n
                hessi(i, j) = (i==j)*(1/(z(i-n)*z(i-n)) + 1. /(precalc(i-n) ^ 2));
            %for p=1:m
            %    hessi(i+m, j+m) = hessi(i+m, j+m) + 1./(precalc(p)*precalc(p));
            %end
            elseif i <= n && j > n
            
            %troisi?me cas : wi zj
            hessi(i,j) = Y(j-n)*X(j-n,i)/(precalc(j-n)^2);
            %hessi(i+m, j) = hessi(i, j+m);
            
            else % i>n et j<=n
                hessi(i,j) = hessi(j,i);
            end
            
            %%% GRAD %%%
            if i <= n && j <= m
                jj = j;
                grad(i) = grad(i) - Y(jj)*X(jj,i)/precalc(jj);
            end
		end
    end
    
    %for i = 1:m
        %X(i,:) = X(i,:)*precalc(i);
    %end
    
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
    while (evalf(X,Y,C,t,w+sz*dx(1:n),z+sz*dx((n+1):((n+m)))) >= evalf(X,Y,C,t,w,z) + alpha*sz*grad'*dx)
        sz=beta*sz;
        %sz
    end
    %x = x+sz*dx;
    w  = w+sz*dx(1:n);
    z  = z+sz*dx((n+1):(n+m));
end

end
