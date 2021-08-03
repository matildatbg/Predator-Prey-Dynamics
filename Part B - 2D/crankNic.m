
%%% The function for the Crank Nicolson scheme, etc. 
%%% Used in the main program for problem B2
function [uh, Mprey] = crankNic(M,A,T,k,p,t)

%Constants 
alpha = 4;
delta = 0.01;
q = delta*k/2;

% Sets inital data for uh
u_initial = zeros(length(p),1);
for j=1:length(p)
    omega = rand();
    u_initial(j,1) = 1+20*omega;
end
xi = u_initial;
uh = xi;

%Sets the initial number of prey
for K = 1:size(t,2)
    nodes = t(1:3,K);
    x = p(1,nodes);
    y = p(2,nodes);
    area_K = polyarea(x,y);
    prey(K) = area_K/3*(xi(nodes(1))+xi(nodes(2))+xi(nodes(3)));
end
Mprey(1) = sum(prey);

%The acutal Crank Nic
time=0;
i = 1;
s = @ (xi) xi./(xi+alpha)-xi.*(1-xi);
while time<=T
    denom = M*xi-k*M*s(xi)-q*A*xi;
    num = M+q*A;
    xi = num\denom;
    uh(:,i+1)= xi;
    
    %Calculates prey
    for K = 1:size(t,2)
        nodes = t(1:3,K);
        x = p(1,nodes);
        y = p(2,nodes);
        area_K = polyarea(x,y);
        prey(K) = area_K/3*(xi(nodes(1))+xi(nodes(2))+xi(nodes(3)));
    end
    Mprey(i+1) = sum(prey);
    i=i+1;
    time=time+k;
end
end
