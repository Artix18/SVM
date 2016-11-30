function [w,z,lambdaDual] = solve(X, Y, C)

m=size(X,1);
n=size(X,2);
epsilon = 0.01;
w=zeros(n,1);
z=2*ones(m,1);
lambdaDual=zeros(n+m,1);

%trouver un pt strict faisable, prendre z tr??s grand ou (w=0 et z>1)

t=1;
mu=10;

while true
	[w2,z2] = central(X,Y,C,t,w,z);
	w=w2;
	z=z2;
    
	if (m/t < epsilon)
		break;
    end
	t = t*mu;
end

for i = 1:n
    lambdaDual(i) = 1/(t*(Y(i)*(X(i,:)*w)+z(i)-1));
end

for i = n+1:(n+m)
    lambdaDual(i) = 1/(t*z(i-n));
end
	
end
