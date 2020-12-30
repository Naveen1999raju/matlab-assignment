T = readtable("../matlab-assignment/data/q3a4.csv");
Top_30 = table2array(T(1:30,2));
x = "Log2FoldChange";
Top_30_genes = table2array(T(1:30,1));
h = heatmap(Top_30);
h.Colormap = autumn;
h.XDisplayLabels = x;
h.YDisplayLabels = Top_30_genes;
ylabel("Genes")