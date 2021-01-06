rng('default') % For reproducibility
n = 200;
x1 = rand(n,1)*2-1;
x2 = rand(n,1)*2-1;
Y = x1-2*x1.*(x2>0)+0.1*rand(n,1);
Mdl = fitrgp([x1 x2],Y);
plotPartialDependence(Mdl,1,'Conditional','centered')