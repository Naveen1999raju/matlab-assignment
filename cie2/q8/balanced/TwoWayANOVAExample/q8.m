load popcorn
popcorn
[p,tbl] = anova2(popcorn,3);
tbl
Fbrands = tbl{2,5}
Fpoppertype = tbl{3,5}
Finteraction = tbl{4,5}