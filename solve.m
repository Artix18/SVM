function [w,z] = solve(X, Y, C)

m=size(X,1);
n=size(X,2);
epsilon = 0.01;
w=zeros(n,1);
z=2*ones(m,1);

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
	
end
