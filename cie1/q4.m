q3a4 = readtable("../matlab-assignment/data/q3a4.csv");
subplot(2,2,1)
quiz4 = q3a4(1:end,:);
scatter(quiz4.B_C_1,quiz4.B_C_2,'r*')
title('scatter plot with 25156 rows, bad presentation')
xlabel('BC1')
ylabel('BC2')

subplot(2,2,2)
quiz4 = q3a4(1:end-10,:);
scatter(quiz4.B_C_1,quiz4.B_C_2,'k*')
title('scatter plot with 25146 rows,its bit ok')
xlabel('BC1')
ylabel('BC2')


subplot(2,2,3)
quiz4 = q3a4(1:end-100,:);
scatter(quiz4.B_C_1,quiz4.B_C_2,'g*')
title('scatter plot with 25056 rows, better than all other')
xlabel('BC1')
ylabel('BC2')


subplot(2,2,4)
quiz4 = q3a4(1:end-1000,:);
scatter(quiz4.B_C_1,quiz4.B_C_2,'b*')
title('scatter plot with 24156 rows, not good')
xlabel('BC1')
ylabel('BC2')
