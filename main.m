m=100; %100 pts
%n=2; %2 dimensions
n=3; %pour avoir un terme affine

%test moyenne al??atoire
%a=1/rand();
%b=1/rand();
%sigma1=[1 0; 0 1];
%sigma2=sigma1;

%test canonique
%a=2; b=2;
%sigma1=[1 0; 0 1];
%sigma2=sigma1;

%test pour l'une des gaussienne ??tal??e vers la premi??re
a=7; b=0.01;
sigma1=[1 0; 0 1];
sigma2=[8 0; 0 1];
%sigma1=sigma2;

pts1=mvnrnd([0;0], sigma1,m/2);
pts2=mvnrnd([a;b], sigma2,m/2);

X = [pts1; pts2];
X = [ones(m,1) X]; %add cst
Y = 2*[zeros(m/2, 1); ones(m/2,1)]-1;

Clist = [0.01; 0.1; 1; 10; 100; 1000; 100000]; %C petit sous-apprend, C grand sur-apprend

for C = Clist'

s = sprintf('C = %f', C);
s
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
S = sprintf('Erreur entrainement : %f \n', (1-nbOk/m)*100);
disp(S);
%nbOk/m*100

figure;
scatter(X(1:m/2,2), X(1:m/2,3), 'o');
hold on;
scatter(X(m/2+1:m,2), X(m/2+1:m,3), 'x');
hold on;
xpts = -20:0.5:20;
droite=-1./w(3) * (w(2)*xpts+w(1));
plot(xpts, droite);
droite2=((a^2+b^2)/2-xpts*a)/b;
plot(xpts, droite2);

title(s);
legend('Donnees classe 1', 'Donnees classe 2', 'Droite du SVM', 'Mediatrice du segment [mu1; mu2]');
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
    if ([1 mvnrnd([0;0], sigma1)]*w<0)
        nbreussi=nbreussi+1;
    end
    if ([1 mvnrnd([a;b], sigma2)]*w>0)
        nbreussi=nbreussi+1;
    end
end
S = sprintf('Erreur de test : %f \n', (1-nbreussi/nbtests)*100);
disp(S);

figure;
semilogy(0:(nbSteps-1), gaps);
title(s);

end
