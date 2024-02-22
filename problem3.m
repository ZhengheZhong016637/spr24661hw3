function problem3
close all;
clear all;

L = 10^4;
y0=sin(pi/4)+10;
Tmax = 10;
figure; hold on;
% maxerr1 = zeros(25,1);
% maxerr2 = zeros(25,1);

% for i=1:25
%     h = 10^(-(1+5*(i-1)/24));
%     Nstep = ceil(Tmax/h);    
%     t=0;
%     u=y0;        
%     maxerr1(i) = log(abs(u-ytrue(t,y0)));
%     for j = 2:Nstep        
%         u = DIRK2(t,u,h);
%         t = t+h;
%         maxerr1(i)=max(maxerr1(i),log(abs(u-ytrue(t,y0))));
%     end
% end
% 
% for i=1:25
%     h = 10^(-(1+5*(i-1)/24));
%     Nstep = ceil(Tmax/h);    
%     t=0;
%     u=y0;        
%     maxerr2(i) = log(abs(u-ytrue(t,y0)));
%     for j = 2:Nstep        
%         u = DIRKo3(t,u,h);
%         t = t+h;
%         maxerr2(i)=max(maxerr2(i),log(abs(u-ytrue(t,y0))));
%     end
% end

% a1 = plot(maxerr1,'-o'); M1 = "DIRK2";
% a2 = plot(maxerr2,'-o'); M2 = "DIRKo3";
% legend([a1,a2],[M1,M2])
% ylabel('log of max error')
% xlabel('h=10^{-(1+5d/24)}')

h=10^-1;
Nstep = ceil(Tmax/h);    
err1 = zeros(Nstep,1);
t1 = 1:Nstep;
t1 = (1/Nstep)*t1;
t=0;
u=y0;        
for j = 2:Nstep        
    u = DIRK2(t,u,h);
    t = t+h;
    err1(j)=log(abs(u-ytrue(t,y0)));
end

h=10^-2;
Nstep = ceil(Tmax/h);    
err2 = zeros(Nstep,1);
t2 = 1:Nstep;
t2 = (1/Nstep)*t2;
t=0;
u=y0;        
for j = 2:Nstep        
    u = DIRK2(t,u,h);
    t = t+h;
    err2(j)=log(abs(u-ytrue(t,y0)));
end

h=10^-3;
Nstep = ceil(Tmax/h);    
err3 = zeros(Nstep,1);
t3 = 1:Nstep;
t3 = (1/Nstep)*t3;
t=0;
u=y0;        
for j = 2:Nstep        
    u = DIRK2(t,u,h);
    t = t+h;
    err3(j)=log(abs(u-ytrue(t,y0)));
end

h=10^-1;
Nstep = ceil(Tmax/h);    
err4 = zeros(Nstep,1);
t4 = 1:Nstep;
t4 = (1/Nstep)*t4;
t=0;
u=y0;        
for j = 2:Nstep        
    u = DIRKo3(t,u,h);
    t = t+h;
    err4(j)=log(abs(u-ytrue(t,y0)));
end

h=10^-2;
Nstep = ceil(Tmax/h);    
err5 = zeros(Nstep,1);
t5 = 1:Nstep;
t5 = (1/Nstep)*t5;
t=0;
u=y0;        
for j = 2:Nstep        
    u = DIRKo3(t,u,h);
    t = t+h;
    err5(j)=log(abs(u-ytrue(t,y0)));
end

h=10^-3;
Nstep = ceil(Tmax/h);    
err6 = zeros(Nstep,1);
t6 = 1:Nstep;
t6 = (1/Nstep)*t6;
t=0;
u=y0;        
for j = 2:Nstep        
    u = DIRKo3(t,u,h);
    t = t+h;
    err6(j)=log(abs(u-ytrue(t,y0)));
end

err1(err1==0) = nan;
err2(err2==0) = nan;
err3(err3==0) = nan;
err4(err4==0) = nan;
err5(err5==0) = nan;
err6(err6==0) = nan;

a1=plot(t1,err1,'-o'); M1="h=10^-1, DIRK2";
a2=plot(t2,err2,'-o'); M2="h=10^-1, DIRK2";
a3=plot(t3,err3,'-o'); M3="h=10^-1, DIRK2";
a4=plot(t4,err4,'-o'); M4="h=10^-1, DIRKo3";
a5=plot(t5,err5,'-o'); M5="h=10^-1, DIRKo3";
a6=plot(t6,err6,'-o'); M6="h=10^-1, DIRKo3";
legend([a1,a2,a3,a4,a5,a6],[M1,M2,M3,M4,M5,M6])
title('T_{max}=10')
ylabel('log of error')
xlabel('time')

    function unp1=DIRK2(t,u,h)
        gamma = 1-1/sqrt(2);
        k1 = (-L*u+L*phi(t+gamma*h)+dphi(t+gamma*h))/(1+L*gamma*h);
        k2 = (-L*u-L*(1-gamma)*h*k1+L*phi(t+h)+dphi(t+h))/(1+L*gamma*h);
        unp1=u+(1-gamma)*h*k1+gamma*h*k2;
    end

    function unp1 = DIRKo3(t,u,h)
        gamma = 1/2+sqrt(3)/6;
        k1 = (-L*u+L*phi(t+gamma*h)+dphi(t+gamma*h))/(1+L*gamma*h);
        k2 = (-L*u-L*(1-2*gamma)*h*k1+L*phi(t+(1-gamma)*h)+dphi(t+(1-gamma)*h))/(1+L*gamma*h);
        unp1 = u+0.5*h*k1+0.5*h*k2;
    end

    function y=ytrue(t,y0)
        y=exp(-L*t)*(y0-phi(0))+phi(t);
    end

    function phi = phi(t)
        phi = sin(t+pi/4);
    end

    function dphi = dphi(t)
        dphi = cos(t+pi/4);
    end
end