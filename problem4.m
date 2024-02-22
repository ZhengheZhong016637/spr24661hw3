[t1,sol1,T1]=RobertsonHW3(1);
[t2,sol2,T2]=RobertsonHW3(2);
[t3,sol3,T3]=RobertsonHW3(3);
% 
% fsz = 20; % fontsize
% figure;hold on
% a1=plot(t1,sol1(:,1),'Linewidth',2); M1 = "h=10^{(-1)}";
% a2=plot(t2,sol2(:,1),'Linewidth',2); M2 = "h=10^{(-2)}";
% a3=plot(t3,sol3(:,1),'Linewidth',2); M3 = "h=10^{(-3)}";
% set(gca,'FontSize',fsz);
% xlabel('t','FontSize',fsz);
% ylabel('x','FontSize',fsz);
% legend([a1,a2,a3],[M1,M2,M3]);
% 
% figure;hold on
% a1=plot(t1,sol1(:,2),'Linewidth',2); M1 = "h=10^{(-1)}";
% a2=plot(t2,sol2(:,2),'Linewidth',2); M2 = "h=10^{(-2)}";
% a3=plot(t3,sol3(:,2),'Linewidth',2); M3 = "h=10^{(-3)}";
% set(gca,'FontSize',fsz);
% xlabel('t','FontSize',fsz);
% ylabel('y','FontSize',fsz);
% legend([a1,a2,a3],[M1,M2,M3]);
% 
% figure;hold on
% a1=plot(t1,sol1(:,3),'Linewidth',2); M1 = "h=10^{(-1)}";
% a2=plot(t2,sol2(:,3),'Linewidth',2); M2 = "h=10^{(-2)}";
% a3=plot(t3,sol3(:,3),'Linewidth',2); M3 = "h=10^{(-3)}";
% set(gca,'FontSize',fsz);
% xlabel('t','FontSize',fsz);
% ylabel('z','FontSize',fsz);
% legend([a1,a2,a3],[M1,M2,M3]);

figure;
bar([T1,T2,T3])