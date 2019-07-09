%Element routine 
function [Ke] = elementroutine(node1, node2)
    DerN = [-1/2, 1/2];
    J = DerN*transpose([node1, node2]);
    Jinv = 1/J;
    B = [-1/2*Jinv, 1/2*Jinv;1/(node1+node2), 1/(node1+node2);1/(node1+node2), 1/(node1+node2)];
    BT = transpose(B);
    C = materialroutine();
    Ke = 2*(BT*C*B)*J*(((node1+node2)/2)^2);
    
end