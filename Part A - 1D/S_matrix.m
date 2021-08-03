
%%% Creates the Stiffness matrix S used in FEM
function S = S_matrix(x,N)
S = zeros(N,N);
for i=1:N-1
    h = x(i+1)-x(i);
    node = [i i+1];
    S(node, node) = S(node, node) + [1 -1; -1 1]/h;
end

% Sets BC to 0
S(1,1) = 10^6; 
S(N,N) = 10^6;
end
