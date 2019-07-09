%material routine
function [C] = materialroutine()
    E = 0.2;                                %Youngs modulus in 
    pratio = 0.20;                          %Poisson ratio
    mu = E/(2*(1+pratio));                  %Mu parameter in given
    lam = pratio*E/((1-2*pratio)*(1+pratio));%Lambda in given
    C = [lam+2*mu lam lam; lam lam+2*mu lam; lam lam lam+2*mu];
end