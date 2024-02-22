clear all;

xmax = 3;
ymax = xmax;
xmin = -xmax;
ymin = -ymax;

x = xmin:0.1:xmax;
y=x;
[X,Y] = meshgrid(x,y);
Z = zeros(size(x,2),size(x,2));
for i = 1:size(x,2)
    for j = 1:size(x,2)
        Z(i,j) = f(x(j),x(i));
    end
end
contourf(X,Y,Z,[1,1],"ShowText",true)
hold on;
plot([0;0],[ymin;ymax])
plot([xmin;xmax],[0;0])
xlabel('Re(h\lambda)')
ylabel('Im(h\lambda)')
axis square;

function f = f(x,y)
    gamma = 1-1/sqrt(2);    
    z = x+1i*y;
    %R = (1+(1-2*gamma)*z)/((1-gamma*z)^2);
    R = (1-2*gamma*z+z)/((1-gamma*z)^2);
    f = abs(R);
end
