load carsmall
X = [Weight,Cylinders,Horsepower];
Y = MPG;
Mdl = fitrtree(X,Y);
view(Mdl,'Mode','graph')
plotPartialDependence(Mdl,1)