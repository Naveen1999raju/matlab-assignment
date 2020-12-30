% read the file 
data = xlsread('../matlab-assignment/data/q5.xlsx');
% remove the last column since it has Nan
data =(data(2:end,1:end-1))';
% assign the labels to value
z = {'1a1m','1a1n','1a1o','1a6z','1a9b','1a9e','1agb','1agc','1agd','1age' ,'1agf' ,'1akj' ,'1ao7' , '1ce6','1cg9','1efx','1exu','1gzp','1hhh','1hsa','1hsb','1i4f','1i3m','1jgd','1jge','6ghn','6ggm' ,'6fgb','6ewo','6ewc','6ewa','6at5','6at9','6bj2','6bj3','6bxp','6c97','6c98','6c99','5vwj','5vz5','5w6a','5w67','5w69','5whk','5wjn','5wkf' ,'5wkh','5wmp'};
% plot the modified data
plot(data)
% add labels title and adjust the legend on the plot 
xlabel('Different sections of gene')
ylabel('Different genes')
title('Relation between different genes(Q 5 )')
legend(z','Location','northeast','NumColumns',10)