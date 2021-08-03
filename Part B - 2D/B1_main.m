
%%% Main program for problem B1
clear all; close all;

geometry = @circleg;
hmax = [1/2 1/4 1/8 1/16 1/32];
EnE = zeros(1,length(hmax));

for i=1:length(hmax)
    [p,e,t] = initmesh(geometry, 'hmax', hmax(i)); %point, edge, element data
    
    u_exact = @(x1,x2) sin(2*pi*x1).*sin(2*pi*x2);
    I = eye(length(p));
    A = stiffness2D(p,t);
    b = loadvect2D(p,t,@RHS);
    
    % BC
    boundaries = zeros(length(e(1,:)),1);
    for j=1:length(e(1,:))
        x1 = p(1,e(1,j));
        x2 = p(2,e(1,j));
        boundaries(j,1) = u_exact(x1,x2);
    end    
    A(e(1 ,:) ,:) = I(e(1 ,:) ,:);
    b(e(1,:),:) = boundaries(:,1);
    
    uh = A\b;
 
    figure(i)
    pdesurf(p,t,uh);
    title("hmax = " + hmax(i));
    
    %calc energy norm of error
    error = u_exact(p(1,:),p(2,:))'-uh;
    EnE(i) = sqrt(error'*A*error); 
end

P = polyfit(log10(hmax), log10(EnE),1);
figure(i+1)
loglog(hmax,EnE,'-b', hmax, hmax.^(P(1)), '-r');
title("Error");
xlabel("log(h_{max})")
ylabel("log(||u-u_h||_E)");
legend("Error", "h_{max}^P, P =" + P(1),'Location', 'southeast');


%Convergence between each hmax
q = zeros(1,length(hmax)-1);
for i=1:length(q)
    q(i) = (log10(EnE(i+1))-log10(EnE(i)))/(log10(hmax(i+1))-log10(hmax(i)));
end
