
%%% Constructs the stiffness matrix
function A = stiffness2D(p,t)

N = size(p,2);
A = sparse(N,N);

for K =1:size(t,2)
    nodes = t(1:3,K);
    x = p(1,nodes);
    y = p(2,nodes);
    area_K = polyarea(x,y);
    bi = [y(2)-y(3); y(3)-y(1); y(1)-y(2)]/(2*area_K);
    ci = [x(3)-x(2); x(1)-x(3); x(2)-x(1)]/(2*area_K);
    AK = (bi*bi'+ci*ci')*area_K;
    A(nodes,nodes) = A(nodes,nodes)+AK;
end

end