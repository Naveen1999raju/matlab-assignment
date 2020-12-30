quiz2 = readtable("../matlab-assignment/data/q2.csv");
boxplot(quiz2.Values,quiz2.Identifiers)
title('Box plot for each type')