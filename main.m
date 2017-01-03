m=100; %100 pts
%n=2; %2 dimensions
n=3; %pour avoir un terme affine
C=1000; %meilleurs resultats avec C grand, (20, 50 ou 100)

a=1/rand();
b=1/rand();
a=2; b=2;
sigma=[1 0; 0 1];

pts1=mvnrnd([a;b], sigma,m/2);
pts2=mvnrnd([0;0], sigma,m/2);

X = [pts1; pts2];
X = [ones(m,1) X]; %add cst
Y = 2*[zeros(m/2, 1); ones(m/2,1)]-1;

[w,z,lambdaDual,Xstep]=solve(X,Y,C);
%[w,z]=solveCVX(X,Y,C);

nbOk = 0;
for i = 1:m
    vv = X(i,:) * w;
   % sg=0;
    if vv >= 0
        sg=1;
    else
        sg=-1;
    end
    if sg==Y(i)
        nbOk = nbOk+1;
    end
end
S = sprintf('Pourcentage correct : %f\n', nbOk/m*100);
disp(S);
%nbOk/m*100

figure(1);
scatter(X(1:m/2,2), X(1:m/2,3), 'o');
hold on;
scatter(X(m/2+1:m,2), X(m/2+1:m,3), 'x');
hold on;
xpts = -20:0.5:20;
droite=-1./w(3) * (w(2)*xpts+w(1));
plot(xpts, droite);
droite2=((a^2+b^2)/2-xpts*a)/b;
plot(xpts, droite2);
hold off;

nbSteps=size(Xstep,1);
gaps = zeros(nbSteps,1);

for i = 1:nbSteps
    wi = Xstep(i, 1:n)';
    zi = Xstep(i, (n+1):(n+m))';
    
    dualityGap = 0.;
    
    for j=1:n
        dualityGap = dualityGap + lambdaDual(j)*(Y(j)*(X(j,:)*wi)+zi(j)-1);
    end
    
    for j=1:m
        dualityGap = dualityGap + lambdaDual(j+n)*zi(j);
    end
    
    gaps(i)=dualityGap;
end

nbtests=10000;
nbreussi=0;
for i=1:nbtests/2
    if ([1 mvnrnd([0;0], sigma)]*w>0)
        nbreussi=nbreussi+1;
    end
    if ([1 mvnrnd([a;b], sigma)]*w<0)
        nbreussi=nbreussi+1;
    end
end
S = sprintf('Pourcentage sur les tests : %f\n', nbreussi*100/nbtests);
disp(S);

figure(2);
plot(0:(nbSteps-1), gaps);
