
%%% Sets the value of f in B1. 
function f = RHS(x1,x2)
    f = 8*pi.^2.*sin(2*pi.*x1).*sin(2*pi.*x2); 
end