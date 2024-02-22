function [t,sol,T]=RobertsonHW3(i)
% Stiff Robertson's problem from chemical kinetics as in
% https://archimede.uniba.it/~testset/report/rober.pdf

% timestep, Tmax, tolearnce for Newton's solver

h = 10^(-i);
Tmax = 1.0e2; % up to 4.0e10
Nsteps = ceil(Tmax/h);
tol = 1.0e-14;
itermax = 20;

y0 = [1.0,0.0,0.0]';

sol = zeros(Nsteps+1,3);
t = h*(1:(Nsteps+1))';
sol(1,:) = y0;

tic % start measuring CPU time

%sol(2,:) = DIRK2step(sol(1,:)',h,tol,itermax)';

for j = 2 : Nsteps
    %sol(j+1,:) = BDF2step(sol(j,:)',sol(j-1,:)',h,tol,itermax)';
    sol(j+1,:)=DIRK2step(sol(j,:)',h,tol,itermax);
end

T=toc;
end
%%

% the right-hand side
function dy = func(y) 
    a = 0.04;
    b = 1.0e4;
    c = 3.0e7;
    dy = zeros(3,1);
    byz = b*y(2)*y(3);
    cy2 = c*y(2)*y(2);
    ax = a*y(1);
    dy(1) = -ax + byz;
    dy(2) = ax - byz - cy2;
    dy(3) = cy2;
end

% the Jacobian matrix for the right-hand side
function J = Jac(y)
    a = 0.04;
    b = 1.0e4;
    c = 3.0e7;
    by = b*y(2);
    bz = b*y(3);
    c2y = 2*c*y(2);
    J = zeros(3);
    J(1,1) = -a;
    J(1,2) = bz;
    J(1,3) = by;
    J(2,1) = a;
    J(2,2) = -bz-c2y;
    J(2,3) = -by;
    J(3,2) = c2y;
end

%% DIRK2

function knew = NewtonIterDIRK2(y,h,k,gamma)
    %note that the newton iteration for DIRK2 and DIRKo3 is the same.
    %since the RHS is time independent we are essentially always solving 
    %k=f(y+gammahk)
    aux = y + h*gamma*k;
    F = k - func(aux);
    DF = eye(3) - h*gamma*Jac(aux);
    knew = k - DF\F;
end

function ynew = NewtonIterBDF2(y,ynm1,ynm2,h)
    F = y - (2/3)*h*func(y)-(4/3)*ynm1+(1/3)*ynm2;
    DF = eye(3)-(2/3)*h*Jac(y);
    ynew = y - DF\F;
end

function ynew = DIRK2step(y,h,tol,itermax)
    gamma = 1.0 - 1.0/sqrt(2);
    k1 = func(y);
    for j = 1 : itermax
        k1 = NewtonIterDIRK2(y,h,k1,gamma);
        if norm(k1 - func(y + h*gamma*k1)) < tol
            break
        end
    end
    k2 = k1;
    y = y + h*(1-gamma)*k1;
    for j =1 : itermax
        k2 = NewtonIterDIRK2(y,h,k2,gamma);
        aux = y + h*gamma*k2;
        if norm(k2 - func(aux)) < tol
            break
        end
    end
    ynew = aux;
end

function ynew = DIRKo3step(y,h,tol,itermax)
    gamma = 1/2+sqrt(3)/6;
    k1 = func(y);
    for j = 1:itermax
        k1 = NewtonIterDIRK2(y,h,k1,gamma);
        if norm(k1-func(y+h*gamma*k1))<tol
            break
        end
    end
    k2=k1;
    aux1=y+h*(1-2*gamma)*k1;
    for j =1 : itermax
        k2 = NewtonIterDIRK2(aux1,h,k2,gamma);
        aux2 = aux1 + h*gamma*k2;
        if norm(k2 - func(aux2)) < tol
            break
        end
    end
    ynew = y+0.5*h*k1+0.5*h*k2;
end

function ynew = BDF2step(y,ynm1,h,tol,itermax)
    ynew = y;
    for j = 1:itermax
        ynew = NewtonIterBDF2(ynew,y,ynm1,h);
        if norm(ynew - (2/3)*h*func(ynew)-(4/3)*y+(1/3)*ynm1)<tol
            break
        end
    end    
end