
%%% Constructs the load vector
function b = loadvect2D(p,t,f)

N = size(p,2);
b = sparse(N,1);

for K = 1:size(t,2)
    nodes = t(1:3,K);
    x = p(1,nodes);
    y = p(2,nodes);
    area_K = polyarea(x,y);
    bK = [f(x(1),y(1)); f(x(2), y(2)); f(x(3),y(3))]/3*area_K;
    b(nodes) = b(nodes)+bK;
end
end