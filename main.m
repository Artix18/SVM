m=100; %100 pts
%n=2; %2 dimensions
C=5; %meilleurs resultats avec C grand, (20, 50 ou 100)

a=1/rand()+100;
b=1/rand()+100;
sigma=[1 0; 0 1];

pts1=mvnrnd([a;b], sigma,m/2);
pts2=mvnrnd([0;0], sigma,m/2);

X = [pts1; pts2];
X = [ones(m,1) X]; %add cst
Y = 2*[zeros(m/2, 1); ones(m/2,1)]-1;

[w,z]=solve(X,Y,C);

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

scatter(X(1:m/2,2), X(1:m/2,3), 'o');
hold on;
scatter(X(m/2+1:m,2), X(m/2+1:m,3), 'x');
hold on;
xpts = -20:0.5:120;
droite=-1./w(3) * (w(2)*xpts+w(1));
plot(xpts, droite);
