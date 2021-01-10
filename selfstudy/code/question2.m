 L=1;
 t=1;
 k=.001;
 n=10;
 nt=500;
 dx=L/n;
 dt=.002;
 alpha=k*dt/dx^2;
 T0=400*ones(1,n);
 T1=300*ones(1,n);
 T0(1) = 300;
 T0(end) = 300;
 for j=1:nt
    for i=2:n-1
        T1(i)=T0(i)+alpha*(T0(i+1)-2*T0(i)+T0(i-1));
    end
    T0=T1;
 end
 dlmwrite('Trial_1.csv',T1,'Delimiter','\n')
 
 figure(1), clf
 plot(T1,'linewidth',2);
 xlabel('Grid points')
 ylabel('Temperature [^oC]')
 title(['Temperature evolution over 10 grid points'])