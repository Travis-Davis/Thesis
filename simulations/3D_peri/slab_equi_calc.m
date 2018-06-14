clear all
clc

%% Parameter Space:
rho_o = 900;
rho_w = 1017;

R_o   = 1;
g     = 9.81;

n_lay = 9;

Vtot  = 4/3*pi*R_o^3;

%% Compute Centroid Locations:
dz = -n_lay*2*R_o;

z = linspace(0,dz,n_lay*100);

for i = 1:length(z)
    for j = 1:n_lay
        pos(i,j)  = z(i) + (2*j-1)*R_o;
    end
end

%% Force Computations:
for i = 1:length(z)
    for j = 1:n_lay
        if     pos(i,j) < -R_o
            Vdisp = Vtot;
        elseif pos(i,j) >= -R_o && pos(i,j) < 0
            Vdisp = Vtot - 1/3*pi*(R_o+pos(i,j))^2*(2*R_o-pos(i,j));
        elseif pos(i,j) >= 0 && pos(i,j) < R_o
            Vdisp = 1/3*pi*(R_o-pos(i,j))^2*(2*R_o+pos(i,j));
        else
            Vdisp = 0;
        end
        
        Fbod(i,j) = rho_w*g*Vdisp - rho_o*g*Vtot;
    end
end

%% Find Contact Force Requirements For Static Equilibrium:
% For static equilibrium, all Fnet must be zero and 3rd Law satisfied:
for i = 1:length(z)
    for j = 1:n_lay
        if     j == 1
            Fcont(i,j) = -Fbod(i,j);
        elseif j == n_lay
            Fcont(i,j) = -Fcont(i,j-1);
        else
            Fcont(i,j) = Fcont(i,j-1)-Fbod(i,j);
        end
    end
end

plot(z,Fcont(:,end),z,Fbod(:,end))
%% Given Contact Force Requirement, find equilibrium position for n'th element:
% Equilibrium occurs when Fbod(i,top) = -Fcont(i,top).
% Find the equilibrium floating element:
diff  = Fbod(:,end)-Fcont(:,end);

z_n   = z(abs(diff) == min(abs(diff)));

n_ele = abs(floor((z_n+R_o)/(2*R_o)));

F_sum = Vtot*g*((n_lay-2*n_ele+2)*rho_o-(n_ele-1)*rho_w);

Vdisp = ((rho_o*n_lay-rho_w*(n_ele-1))/rho_w)*Vtot;

F_c_below = (n_ele-1)*Vtot*(rho_w-rho_o)*g;
F_c_above = (n_lay-n_ele)*Vtot*rho_o*g;
F_b_sing  = Vtot*(rho_w-rho_o)*g;


%% Test:
n_test = find(ones(1,n_lay));
Vdisp = (rho_o/rho_w*n_lay-n_test)*Vtot;

figure(1)
plot([0,n_test(end)],[0,0])
hold on
plot([0,n_test(end)],[Vtot,Vtot])
scatter(n_test,Vdisp)