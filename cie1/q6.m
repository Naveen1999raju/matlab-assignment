T = readtable('../matlab-assignment/data/q6.xlsx');
figure
T.ATC=abs(fft(T.CH00));
plot(T.f,T.ATC)
title("frequency vs Atc")
T.ATC1=abs(fft(T.CH01));
figure
plot(T.f,T.ATC1)
title("frequency vs Atc1")
T.ATC2=abs(fft(T.CH02));
figure
plot(T.f,T.ATC2)
title("frequency vs Atc2")
figure
plot(T.f,T.ATC,T.ATC1,T.ATC2)
title("frequecy ,ATC,ATC1,ATC2")