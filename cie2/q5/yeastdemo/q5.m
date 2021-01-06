load yeastdata.mat
numel(genes)
genes{15}

plot(times, yeastvalues(15,:))
xlabel('Time (Hours)');
ylabel('Log2 Relative Expression Level');

plot(times, 2.^yeastvalues(15,:))
xlabel('Time (Hours)');
ylabel('Relative Expression Level');

hold on
plot(times, 2.^yeastvalues(16:26,:)')
xlabel('Time (Hours)');
ylabel('Relative Expression Level');
title('Profile Expression Levels');

emptySpots = strcmp('EMPTY',genes);
yeastvalues(emptySpots,:) = [];
genes(emptySpots) = [];
numel(genes)

nanIndices = any(isnan(yeastvalues),2);
yeastvalues(nanIndices,:) = [];
genes(nanIndices) = [];
numel(genes)

mask = genevarfilter(yeastvalues);
 
yeastvalues = yeastvalues(mask,:);
genes = genes(mask);
numel(genes)


[mask,yeastvalues,genes] = genelowvalfilter(yeastvalues,genes,'absval',log2(3));
numel(genes)

[mask,yeastvalues,genes] = geneentropyfilter(yeastvalues,genes,'prctile',15);
numel(genes)
corrDist = pdist(yeastvalues,'corr');
clusterTree = linkage(corrDist,'average');

clusters = cluster(clusterTree,'maxclust',16);

figure
for c = 1:16
    subplot(4,4,c);
    plot(times,yeastvalues((clusters == c),:)');
    axis tight
end
suptitle('Hierarchical Clustering of Profiles');

rng('default');

[cidx, ctrs] = kmeans(yeastvalues,16,'dist','corr','rep',5,'disp','final');
figure
for c = 1:16
    subplot(4,4,c);
    plot(times,yeastvalues((cidx == c),:)');
    axis tight
end
suptitle('K-Means Clustering of Profiles');


figure
for c = 1:16
    subplot(4,4,c);
    plot(times,ctrs(c,:)');
    axis tight
    axis off   
end
suptitle('K-Means Clustering of Profiles');

cgObj = clustergram(yeastvalues(:,2:end),'RowLabels',genes,'ColumnLabels',times(2:end));
h = mapcaplot(yeastvalues,genes);

[pc, zscores, pcvars] = pca(yeastvalues);

pcvars./sum(pcvars) * 100

cumsum(pcvars./sum(pcvars) * 100)

figure
scatter(zscores(:,1),zscores(:,2));
xlabel('First Principal Component');
ylabel('Second Principal Component');
title('Principal Component Scatter Plot');

figure
pcclusters = clusterdata(zscores(:,1:2),'maxclust',8,'linkage','av');
gscatter(zscores(:,1),zscores(:,2),pcclusters)
xlabel('First Principal Component');
ylabel('Second Principal Component');
title('Principal Component Scatter Plot with Colored Clusters');

if ~exist(which('selforgmap'),'file')
    disp('The Self-Organizing Maps section of this example requires the Deep Learning Toolbox.')
    return
end

P = zscores(:,1:2)';
net = selforgmap([4 4]);
net = train(net,P);
figure
plot(P(1,:),P(2,:),'.g','markersize',20)
hold on
plotsom(net.iw{1,1},net.layers{1}.distances)
hold off

distances = dist(P',net.IW{1}');
[d,cndx] = min(distances,[],2); % cndx contains the cluster index

figure
gscatter(P(1,:),P(2,:),cndx); legend off;
hold on
plotsom(net.iw{1,1},net.layers{1}.distances);
hold off

close('all');
delete(cgObj);
delete(h);
nntraintool('close');
