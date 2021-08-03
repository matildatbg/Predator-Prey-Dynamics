

%%% Sets the RHS of the second degree ODE. 
function f = RHSfunc(x)
f = zeros(1, length(x));
    for i=1:length(x)
        if abs(x(i)) <= 0.1
            f(i) = -5;
        elseif abs(0.2-abs(x(i))) <=0.1
            f(i) = 25;
        elseif abs(0.6-abs(x(i))) <=0.1
            f(i) = -30;
        elseif abs(0.9-abs(x(i))) <=0.2
            f(i) = 20;
        else
            f(i) = 1;
        end
    end
end