m=100; %100 pts
n=2; %2 dimensions
C=0.01;

a=1/rand()+10;
b=1/rand()+10;
sigma=[1 0; 0 1];

pts1=mvnrnd([a;b], sigma,m/2);
pts2=mvnrnd([0;0], sigma,m/2);

X = [pts1; pts2];
Y = floor(rand(m,1)+0.5)*2-1;
%Y(2)=-Y(1);

w=solve(X,Y,C);

nbOk = 0;
for i = 1:m
    vv = X(i,:) * w;
    sg=0;
    if vv >= 0
        sg=1;
    else
        sg=-1;
    end
    if sg==Y(i)
        nbOk = nbOk+1;
    end
end
nbOk
nbOk/m*100
