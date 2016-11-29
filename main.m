m=100; %100 pts
n=2; %2 dimensions
C=5; %meilleurs resultats avec C grand, (20, 50 ou 100)

a=1/rand();
b=1/rand();
sigma=[1 0; 0 1];

pts1=mvnrnd([a;b], sigma,m/2);
pts2=mvnrnd([0;0], sigma,m/2);

X = [pts1; pts2];
Y = 2*[zeros(m/2, 1); ones(m/2,1)]-1;

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
