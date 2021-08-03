

%%% The main program for problem B2. 
%%% Change alpha and delta in crankNic function
clear all; close all;

% Geometry setup
geometry = @circleg;
hmax = [1/5 1/20 1/40];

% Constants etc
T = 2;
k = 0.01*hmax;
uh = cell(1,length(hmax));
Mprey = cell(1,length(hmax));
p = cell(1,length(hmax));
t = cell(1,length(hmax));

for i=1:length(hmax)
    [p{1,i},e,t{1,i}] = initmesh(geometry, 'hmax', hmax(i));
    A = stiffness2D(p{1,i},t{1,i});
    M = mass2D(p{1,i},t{1,i});
    
    %Crank-Nicolson
    [uh{1,i}, Mprey{1,i}] = crankNic(M,A,T,k(i),p{1,i},t{1,i});
    
    %Plots uh
    figure(i)
    subplot(1,2,1)
    pdesurf(p{1,i}, t{1,i}, uh{1,i}(:,1));
    title("T = 0")
    subplot(1,2,2)
    pdesurf(p{1,i}, t{1,i}, uh{1,i}(:,end));
    title("T = " + T);
    suptitle("h_{max} = " + hmax(i))
end

%Plots Mprey
figure(i+1)
for j=1:length(hmax)
    time = 0:k(j):T;
    subplot(3,1,j);
    plot(time, Mprey{1,j})
    title("h_{max} = " + hmax(j));
    legend("Mprey");
    suptitle("Numb. of prey over time")
end

