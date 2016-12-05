function [w,z] = solveCVX(X,Y,C)
%on utilise ici CVX
m=size(X,1);
n=size(X,2);

cvx_begin
    variable w(n)
    variable z(m)
    minimize( 0.5*w'*w + C*ones(1,m)*z )
    subject to
        Y.*(X*w) >= ones(m,1)-z
        z >= zeros(m,1)
cvx_end

end

