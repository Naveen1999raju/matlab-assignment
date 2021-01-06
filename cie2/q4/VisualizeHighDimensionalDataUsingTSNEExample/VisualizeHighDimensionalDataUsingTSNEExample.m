%% Visualize High-Dimensional Data Using t-SNE
% This example shows how to visualize the MNIST data [1], which consists of
% images of handwritten digits, using the |tsne| function. The images are
% 28-by-28 pixels in grayscale. Each image has an associated label from 0
% through 9, which is the digit that the image represents. |tsne| reduces
% the dimension of the data from 784 original dimensions to 50 using PCA,
% and then to two or three using the t-SNE Barnes-Hut algorithm.
%% Obtain Data
% Begin by obtaining image and label data from
%
% http://yann.lecun.com/exdb/mnist/
%
% Unzip the files. For this example, use the |t10k-images| data.
imageFileName = 't10k-images-idx3-ubyte/s1036Umf9eAy';
labelFileName = 't10k-labels-idx1-ubyte/s103mWIWvWwn';
%%
% Process the files to load them in the workspace. The code for this
% processing function appears at the end of this example.
[X,L] = processMNISTdata(imageFileName,labelFileName);
%% Reduce Dimension of Data to Two
% Obtain two-dimensional analogues of the data clusters using t-SNE. Use
% PCA to reduce the initial dimensionality to 50. Use the Barnes-Hut
% variant of the t-SNE algorithm to save time on this relatively large data
% set.
rng default % for reproducibility
Y = tsne(X,'Algorithm','barneshut','NumPCAComponents',50);
%%
% Display the result, colored with the correct labels.
figure
gscatter(Y(:,1),Y(:,2),L)
%%
% t-SNE creates clusters of points based solely on their relative
% similarities that correspond closely to the true labels.
%% Reduce Dimension of Data to Three
% t-SNE can also reduce the data to three dimensions. Set the |tsne|
% |'NumDimensions'| name-value pair to |3|.
rng default % for fair comparison
Y3 = tsne(X,'Algorithm','barneshut','NumPCAComponents',50,'NumDimensions',3);
figure
scatter3(Y3(:,1),Y3(:,2),Y3(:,3),15,L,'filled');
view(-93,14)
%%
% Here is the code of the function that reads the data into the workspace.
%
% <include>processMNISTdata.m</include>
%
%% References
% [1] Yann LeCun (Courant Institute, NYU) and Corinna Cortes (Google Labs,
% New York) hold the copyright of MNIST dataset, which is a derivative work
% from original NIST datasets. MNIST dataset is made available under the
% terms of the Creative Commons Attribution-Share Alike 3.0 license,
% https://creativecommons.org/licenses/by-sa/3.0/
%% 
% Copyright 2012-2019 The MathWorks, Inc.