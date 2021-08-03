


%%% Creates the Mass matrix M used in FEM
function M = M_matrix(x,N)
M = zeros(N,N);

for i=1:N-1
    h = x(i+1)-x(i);
    node = [i i+1];
    M(node,node) = M(node, node) + [2 1; 1 2]*h/6;
end
end