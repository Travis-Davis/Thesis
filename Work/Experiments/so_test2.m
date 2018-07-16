clear all
clc

Go = .4;
K = 5.81E9;

tol = .01;
hor = 1;
L = linspace(1/10*hor-tol,hor-tol,1001);

for P = 1:length(L)
    x = -(ceil(hor/L(P))+1):L(P):(ceil(hor/L(P))+1);
    y = x;
    z = x;
    
    ind = 0;
    xi_1 = [];
    while ind*L(P) > L(P)-hor
        for i = 1:length(x)
            for j = 1:length(y)
                for k = 1:length(z)
                    xi = sqrt(x(i)^2+y(j)^2+(z(k)-ind))^2;
                    
                    if xi < hor && z(k) > 0
                        xi_1(end+1) = xi;
                    end
                end
            end
        end
        ind = ind-1;
    end
    
    c = 18*K/(pi*hor^4);
    s0_1(P) = sqrt((2*Go)/(c*sum(xi_1.^2)));
    disp(L(P)) 
end

s00 = sqrt((5*Go)/(9*K*hor));
%%
plot(L,s00)
hold on
scatter(L,s0_1)
refline(0,s00)
xlim([0,max(hor)-tol])
ylim([0,max([max(s0_1),max(s00)])])