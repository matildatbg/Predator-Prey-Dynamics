

%%% Main script for part A.
%%% Finds the solution for u_h by using FEM, caluclates
%%% a posteriori error estimation, and then runs the loop
%%& until desired error tolerance is reached.

close all; clear all;
a = -1;
b = 1;
delta = 0.01;
N = 12;
x = linspace(a,b,N);
stopTOL = 10^-3;
stopNodes = 10^4;
condition = 1;

while (condition > stopTOL && N < stopNodes)
    %FEM
    f = RHSfunc(x)./(delta);
    bLoad = load_vect(x,f);
    S = S_matrix(x,length(x));
    M = M_matrix(x,length(x));
    uh = (S\bLoad);         
    uhLap = -(M\(S*uh));
    
    %Calculates residual and eta^2
    Rh = delta.*f+delta.*uhLap';
    eta2 = zeros(N-1,1);
    for i=1:N-1
        h = x(i+1)-x(i);
        temp = h/2*(Rh(i)^2+Rh(i+1)^2);
        eta2(i) = h^2*temp;
    end
    
    %Refines the mesh
    lambda = 0.9;
    xold = x;
    for i = 1:length(eta2)
        if eta2(i) > lambda*max(eta2)
            x = [x (x(i+1)+x(i))/2];
        end
    end
    x = sort(x);
    N = length(x);
    %condition = sum(eta2);
    condition = sum(sqrt(eta2));
end

% Plots the results
subplot(2,2,1)
plot(xold,uh);
title("Solution u_h");
xlabel("x");

subplot(2,2,2)
plot(xold,Rh);
title("Residual");
xlabel("x");
ylabel("R(u_h)");

subplot(2,2,3)
plot(xold(2:end),sqrt(eta2));
title("Error");
xlabel("x");
ylabel("\eta(u_h)");

subplot(2,2,4)
plot(xold(2:end),[1./diff(xold)])
title("Mesh Distribution");
xlabel("x");
ylabel("Amount");

