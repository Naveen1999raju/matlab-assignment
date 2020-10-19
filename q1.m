quiz1 = readtable("../matlab-assignment/data/q1.csv");
width = 0.7;
lim = 500;
figure
set(gcf, 'numbertitle', 'off', 'name', 'radialBar demo')
ax = subplot(1,1,1,'box', 'on');
radialBar(ax, quiz1.NumberOfPeptides, width, [pi/2 -2*pi]);
title("Radialbar of Number of Peptides of different species")
labels = {'Atg8 protein family ligands (4)','Phosphrylation (27)','SUMO interaction site (3)','SH3 motif (5)','Ubiquitination/deubiquitinating enzyme (4)','SH2 domains binding motif (3)','PTM of histones (4)','N-glycosylation site (3)','IAP-binding motif (3)'};
figure
pie(quiz1.NumberOfPeptides)
legend(labels,'Location','southeast')
