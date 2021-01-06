load popcorn
popcorn
[~,~,stats] = anova2(popcorn,3,'off')
c = multcompare(stats)
c = multcompare(stats,'Estimate','row')