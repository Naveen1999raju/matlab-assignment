assignment = xlsread("../matlab-assignment/data/assignment.xlsx");
assign1 = assignment(1 , 1:96);
assign2 = assignment(2 , 1:96);
hold on
plot(assign1,'b','linewidth',2);
plot(assign2,'r','linewidth',2);
hold off
xlabel("Different mutation sites")
legend('sum','occur','FontSize',12)
title("Mutation")