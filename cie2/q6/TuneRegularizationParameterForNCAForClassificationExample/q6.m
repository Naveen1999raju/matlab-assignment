load('twodimclassdata.mat')
figure
gscatter(X(:,1),X(:,2),y)
xlabel('x1')
ylabel('x2')
n = size(X,1);
rng('default')
XwithBadFeatures = [X,randn(n,100)*sqrt(20)];
XwithBadFeatures = (XwithBadFeatures-min(XwithBadFeatures,[],1))./range(XwithBadFeatures,1);
X = XwithBadFeatures;
ncaMdl = fscnca(X,y,'FitMethod','exact','Verbose',1, ...
    'Solver','lbfgs');
figure
semilogx(ncaMdl.FeatureWeights,'ro')
xlabel('Feature index')
ylabel('Feature weight')   
grid on
cvp           = cvpartition(y,'kfold',5);
numtestsets   = cvp.NumTestSets;
lambdavalues  = linspace(0,2,20)/length(y); 
lossvalues    = zeros(length(lambdavalues),numtestsets);
for i = 1:length(lambdavalues)                
    for k = 1:numtestsets
        
        % Extract the training set from the partition object
        Xtrain = X(cvp.training(k),:);
        ytrain = y(cvp.training(k),:);
        
        % Extract the test set from the partition object
        Xtest = X(cvp.test(k),:);
        ytest = y(cvp.test(k),:);
        
        % Train an nca model for classification using the training set
        ncaMdl = fscnca(Xtrain,ytrain,'FitMethod','exact', ...
            'Solver','lbfgs','Lambda',lambdavalues(i));
        
        % Compute the classification loss for the test set using the nca
        % model
        lossvalues(i,k) = loss(ncaMdl,Xtest,ytest, ...
            'LossFunction','quadratic');   
   
    end                          
end
figure
plot(lambdavalues,mean(lossvalues,2),'ro-')
xlabel('Lambda values')
ylabel('Loss values')
grid on
[~,idx] = min(mean(lossvalues,2)); % Find the index
bestlambda = lambdavalues(idx) % Find the best lambda value
ncaMdl = fscnca(X,y,'FitMethod','exact','Verbose',1, ...
     'Solver','lbfgs','Lambda',bestlambda);
 figure
semilogx(ncaMdl.FeatureWeights,'ro')
xlabel('Feature index')
ylabel('Feature weight')    
grid on