clear all
clc

Go = .4;
K = 5.81E9;

tol = .01;
L = 3/5;
hor = linspace(L+tol,20+tol,1001);

for P = 1:length(hor)
    x = -(ceil(hor(P)/L)+1):L:(ceil(hor(P)/L)+1);
    y = x;
    z = x;
    
    ind = 0;
    xi_1 = [];
    while ind > L-hor(P)
        for i = 1:length(x)
            for j = 1:length(y)
                for k = 1:length(z)
                    xi = sqrt(x(i)^2+y(j)^2+(z(k)-ind))^2;
                    
                    if xi < hor(P) && z(k) > 0
                        xi_1(end+1) = xi;
                    end
                end
            end
        end
        ind = ind-1;
    end
    
    c = 18*K/(pi*hor(P)^4);
    s0_1(P) = sqrt((2*Go)/(c*sum(xi_1.^2)*L^2));
    s00(P) = sqrt((5*Go)/(9*K*hor(P)));
    disp(hor(P))
end

%%
plot(hor,s00)
hold on
scatter(hor,s0_1)
plot([5,5],[0,max(s00)],'--k','LineWidth',3)
set(gca,'FontSize',24)
xlim([0,max(hor)-tol])
xlabel('$\delta$ (m)')
ylabel('Critical Strain, s00')
title({'$\delta$ vs. s00','($L = 0.6 m, G_{o} = 0.4 J/m^{2}, K = 5.81 GPa$)'}) 
ylim([0,max(s00)])
legend('Continuous Solid','Discrete: CP=(0,0,1)','Neighbor Overflow Limit')