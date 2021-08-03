
%%% Creates the load vector used in FEM
%%% Code obtained from instructions to Lab1 in the
%%% course "Applied Finite Element methods" @ Uppsala University
function B=load_vect(x,RHS)
N = length(x);
B = zeros(N, 1);
for i = 1:N-1
    h = x(i+1) - x(i);
    n = [i i+1];
    B(n) = B(n) + [RHS(i); RHS(i+1)]*h/2;
end