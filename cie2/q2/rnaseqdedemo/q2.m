
load pasilla_count_noMM.mat

geneCountTable(1:10,:)


exonCountTable(1:10,:)
samples = geneCountTable(:,3:end).Properties.VariableNames;
untreated = strncmp(samples,'untreated',length('untreated'))
treated = strncmp(samples,'treated',length('treated'))
st = geneSummaryTable({'Assigned','Unassigned_ambiguous','Unassigned_noFeature'},:)
bar(table2array(st)','stacked');
legend(st.Properties.RowNames','Interpreter','none','Location','southeast');
xlabel('Sample')
ylabel('Number of reads')
geneSummaryTable{'TotalEntries',:} - sum(geneSummaryTable{2:end,:})
chr2L = strcmp(exonCountTable.Reference,'2L');
exonCount = exonCountTable{:,3:end};
exonStart = regexp(exonCountTable{chr2L,1},'_(\d+)_','tokens');
exonStart = [exonStart{:}];
exonStart = cellfun(@str2num, [exonStart{:}]'); 
[~,idx] = sort(exonStart);

% plot read coverage along the genomic coordinates
figure;
plot(exonStart(idx),exonCount(idx,treated),'.-r',...
exonStart(idx),exonCount(idx,untreated),'.-b');
xlabel('Genomic position');
ylabel('Read count (exon level)');
title('Read count on Chromosome arm 2L');

% plot read coverage for each group separately 
figure;
subplot(2,1,1); 
plot(exonStart(idx),exonCount(idx,untreated),'.-r');
ylabel('Read count (exon level)');
title('Chromosome arm 2L (treated samples)');
subplot(2,1,2); 
plot(exonStart(idx),exonCount(idx,treated),'.-b');
ylabel('Read count (exon level)');
xlabel('Genomic position');
title('Chromosome arm 2L (untreated samples)');

load pasilla_geneLength
geneLength(1:10,:)

chr3 = (geneLength.Reference == '3L') | (geneLength.Reference == '3R');
sum(chr3) 
counts = geneCountTable{:,3:end};
[~,j,k] = intersect(geneCountTable{:,'ID'},geneLength{chr3,'ID'});
gstart = geneLength{k,'Start'};
gcounts = counts(j,:);

% sort according to ascending start position
[~,idx] = sort(gstart);

% plot read coverage by genomic position
figure;
plot(gstart(idx), gcounts(idx,treated),'.-r',...
	gstart(idx), gcounts(idx,untreated),'.-b');
xlabel('Genomic position')
ylabel('Read count');
title('Read count on Chromosome 3');
% estimate pseudo-reference with geometric mean row by row
pseudoRefSample = geomean(counts,2); 
nz = pseudoRefSample > 0; 
ratios = bsxfun(@rdivide,counts(nz,:),pseudoRefSample(nz));
sizeFactors = median(ratios,1)

% transform to common scale
normCounts = bsxfun(@rdivide,counts,sizeFactors);
normCounts(1:10,:)

figure;

subplot(2,1,1)
maboxplot(log2(counts),'title','Raw Read Count','orientation','horizontal')
ylabel('sample')
xlabel('log2(counts)')

subplot(2,1,2)
maboxplot(log2(normCounts),'title','Normalized Read Count','orientation','horizontal')
ylabel('sample')
xlabel('log2(counts)')


% consider the mean
meanTreated = mean(normCounts(:,treated),2); 
meanUntreated = mean(normCounts(:,untreated),2);

% consider the dispersion
dispTreated = std(normCounts(:,treated),0,2) ./ meanTreated; 
dispUntreated = std(normCounts(:,untreated),0,2) ./ meanUntreated;

% plot on a log-log scale
figure;
loglog(meanTreated,dispTreated,'r.');
hold on;
loglog(meanUntreated,dispUntreated,'b.');
xlabel('log2(Mean)');
ylabel('log2(Dispersion)');
legend('Treated','Untreated','Location','southwest');


% compute the mean and the log2FC
meanBase = (meanTreated + meanUntreated) / 2;
foldChange = meanTreated ./ meanUntreated;
log2FC = log2(foldChange);

% plot mean vs. fold change (MA plot)
mairplot(meanTreated, meanUntreated,'Type','MA','Plotonly',true);
set(get(gca,'Xlabel'),'String','mean of normalized counts')
set(get(gca,'Ylabel'),'String','log2(fold change)')


mairplot(meanTreated,meanUntreated,'Labels',geneCountTable.ID,'Type','MA');

geneTable = table(meanBase,meanTreated,meanUntreated,foldChange,log2FC);
geneTable.Properties.RowNames = geneCountTable.ID;

%%

% summary 
summary(geneTable)

%%

% preview 
geneTable(1:10,:)

%%

% access information about a specific gene
myGene = 'FBgn0261570';
geneTable(myGene,:)
geneTable(myGene,{'meanBase','log2FC'})

% access information about a given gene list
myGeneSet = {'FBgn0261570','FBgn0261573','FBgn0261575','FBgn0261560'};
geneTable(myGeneSet,:)

tIdentity = nbintest(counts(:,treated),counts(:,untreated),'VarianceLink','Identity');
h = plotVarianceLink(tIdentity);

% set custom title
h(1).Title.String = 'Variance Link on Treated Samples';
h(2).Title.String = 'Variance Link on Untreated Samples';

tConstant = nbintest(counts(:,treated),counts(:,untreated),'VarianceLink','Constant');
h = plotVarianceLink(tConstant);

% set custom title
h(1).Title.String = 'Variance Link on Treated Samples';
h(2).Title.String = 'Variance Link on Untreated Samples';

tLocal = nbintest(counts(:,treated),counts(:,untreated),'VarianceLink','LocalRegression');
h = plotVarianceLink(tLocal);

% set custom title
h(1).Title.String = 'Variance Link on Treated Samples';
h(2).Title.String = 'Variance Link on Untreated Samples';

h = plotVarianceLink(tLocal,'compare',true);

% set custom title
h(1).Title.String = 'Variance Link on Treated Samples';
h(2).Title.String = 'Variance Link on Untreated Samples';
  

figure;
histogram(tLocal.pValue,100)
xlabel('P-value')
ylabel('Frequency')
title('P-value enrichment')


lowCountThreshold = 10;
lowCountGenes = all(counts < lowCountThreshold, 2);
histogram(tLocal.pValue(~lowCountGenes),100)
xlabel('P-value')
ylabel('Frequency')
title('P-value enrichment without low count genes')



% compute the adjusted P-values (BH correction)
padj = mafdr(tLocal.pValue,'BHFDR',true);

% add to the existing table
geneTable.pvalue = tLocal.pValue;
geneTable.padj = padj;

% create a table with significant genes
sig = geneTable.padj < 0.1;
geneTableSig = geneTable(sig,:);
geneTableSig = sortrows(geneTableSig,'padj');
numberSigGenes = size(geneTableSig,1)

% find up-regulated genes
up = geneTableSig.log2FC > 1;
upGenes = sortrows(geneTableSig(up,:),'log2FC','descend');
numberSigGenesUp = sum(up)

% display the top 10 up-regulated genes
top10GenesUp = upGenes(1:10,:)

% find down-regulated genes
down = geneTableSig.log2FC < -1;
downGenes = sortrows(geneTableSig(down,:),'log2FC','ascend');
numberSigGenesDown = sum(down)

% find top 10 down-regulated genes
top10GenesDown = downGenes(1:10,:)



figure
scatter(log2(geneTable.meanBase),geneTable.log2FC,3,geneTable.padj,'o')
colormap(flipud(cool(256)))
colorbar;
ylabel('log2(Fold Change)')
xlabel('log2(Mean of normalized counts)')
title('Fold change by FDR')

