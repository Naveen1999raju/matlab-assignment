quiz9 = readtable("../matlab-assignment/data/q9.xlsx");
T_array = table2array(quiz9);
xnames = quiz9.Properties.VariableNames;
gplotmatrix(T_array,[],[], [],[],[],[],'variable',xnames);