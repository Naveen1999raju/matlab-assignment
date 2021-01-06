imageFileName = 't10k-images-idx3-ubyte/s1036Umf9eAy';
labelFileName = 't10k-labels-idx1-ubyte/s103mWIWvWwn';

[X,L] = processMNISTdata(imageFileName,labelFileName);

rng default % for reproducibility
Y = tsne(X,'Algorithm','barneshut','NumPCAComponents',50);
figure
gscatter(Y(:,1),Y(:,2),L)

Y3 = tsne(X,'Algorithm','barneshut','NumPCAComponents',50,'NumDimensions',3);
figure
scatter3(Y3(:,1),Y3(:,2),Y3(:,3),15,L,'filled');
view(-93,14)
