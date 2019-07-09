%Assignment

% Given Data

tau = 0.01;             %time like parameter
Initial_stress = 200;   %MPa
ir = 5;                 %inner radius
or = 20;                %outer radius
Epsilon = 0;            %Guass point
del_u = 1;              %Error
n = 50;                 %no of elements
volstr = -0.005;        %volumetric strain


%Mesh Generation

meshrefinementfactor=5; %ratio of element sizes at outer and inner radius
q=meshrefinementfactor^(1./(n-1));%ratio between element sizes of subsequent elements for a geometric series
lelem=(or-ir)*(1-q)/(1-meshrefinementfactor*q);%size of first element
rnode=ir;
rnodes=[ir];
%loop over all elements
for i=1:n
        rnode=rnode+lelem;
        rnodes=[rnodes;rnode];
        lelem=lelem*q;
end

%Initial displacement guesses

glo_disp = transpose(linspace((1/3*tau*(-volstr)*rnodes(1)), 0, n+1));
red_disp = glo_disp(2:end);

%Main Program
while norm(del_u)>(0.005*norm(red_disp))
    glo_k = zeros(n+1);
    for i=1:n
        Ae = zeros(2,n+1);
        Ae(1,i)=1;
        Ae(2,i+1)=1;
        AeT = transpose(Ae);
        Ke = elementroutine(rnodes(i),rnodes(i+1));
        Ke = AeT*Ke*Ae;
        glo_k = glo_k + Ke;
    end

    kred = glo_k(2:end, 2:end);
    glo_fext = zeros(n+1, 1);

    %Newton Raphson method

    G = glo_k*glo_disp-glo_fext;
    Gred = G(2:end);
    del_u = kred\Gred;  %[Dis,r]=linsolve(red_Kg,red_F_ex) %different way of solving
    red_disp = red_disp - del_u;

    ini_disp = (1/3*tau*(-volstr)*rnodes(1));
    glo_disp = [ini_disp; red_disp];
end
   
%actual elastic strain
u_exact = ((ir^3*(-volstr)*tau)./(3*(rnodes.^2)));

%Plotting 

figure
plot(rnodes, glo_disp, rnodes, u_exact, '--')
title('Exact Vs Analytical Displacement')
xlabel('Nodal coordinates')
ylabel('Displacement')
legend('Analytical', 'Exact')
