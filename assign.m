assignment = xlsread("../matlab-assignment/data/statement.xlsx");
assignment1 = assignment(139,1:96);

plot(assignment1,'linewidth',2);
xlabel("Different mutation sites")
title("Mutation")
