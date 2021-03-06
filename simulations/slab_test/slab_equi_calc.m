clear all
clc

platform = 'mb';
input    = 'slab_test';
%% Parameter Space:
x_lo = -50;
x_hi = 50;
y_lo = -50;
y_hi = 50;

rho_o = 900;
rho_w = 1017;

R_o   = 1;
g     = 9.81;

n_lay = 1;

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

%% Given Contact Force Requirement, find equilibrium position for n'th element:
% Equilibrium occurs when Fbod(i,top) = -Fcont(i,top).

if n_lay == 1;
    cub = roots([rho_w*pi/3,rho_w*2*pi*R_o,-rho_w*pi*R_o^2,-rho_o*4/3*pi*R_o^3]);
    for i = 1:3
        if cub(i) < 0
            z_bot = cub(i) - R_o;
            z_top = z_bot + 2*R_o;
        end
    end
else
    
    n_test = find(ones(1,n_lay));
    Vdisp  = (rho_o/rho_w*n_lay-n_test)*Vtot;

    n_b = find(Vdisp <= Vtot & Vdisp > 0)-1;
    n_a = n_lay - n_b - 1;

    figure
    scatter(n_test,Vdisp)
    refline([0 0])
    refline([0 Vtot])

    if Vdisp(n_b+1) <= 1/2*Vtot
        Vfoo = Vdisp(n_b+1);
    else
        Vfoo = Vtot - Vdisp(n_b+1);
    end

    cub = roots([1,-3*R_o,0,3/pi*Vfoo]);
    for i = 1:3
        if cub(i) <= R_o && cub(i) > 0 
            if Vdisp(n_b+1) <= 1/2*Vtot
                h = cub(i);
                z_bot = -(2*R_o*n_b + h)-R_o;
            else
                h = cub(i);
                z_bot = -((2*R_o-h) + 2*R_o*n_b)-R_o;
            end
            z_top = n_lay*2*R_o+z_bot;
        end
    end
end

%% Run the Simulation:
% Set required Parameters:
m = Vtot*rho_o;
n_x = floor((x_hi - x_lo)/(2*R_o));
n_y = floor((y_hi - y_lo)/(2*R_o));

system('rm atom.include *.log log.* ../data/DUMP/*.vtk ../data/DUMP/*.dump');

fileID = fopen('atom.include','w');
for i = 1:n_x
    for j = 1:n_y
        for k = 1:n_lay
            fprintf(fileID,sprintf('create_atoms 1 single %0.32f %0.32f %0.32f\n',...
                x_lo+(i-1)*2*R_o,y_lo+(j-1)*2*R_o,z_bot+(k-1)*2*R_o));
        end
    end
end

fprintf(fileID,sprintf('set group all volume %0.32f\n',Vtot));
fprintf(fileID,'set group all density  ${dens}\n');
fprintf(fileID,sprintf('set group all mass %0.32f\n',m));

fprintf(fileID,sprintf('variable diam equal %0.32f\n',2*R_o));
fprintf(fileID,sprintf('variable rad equal ${diam}*0.5\n'));
fprintf(fileID,sprintf('variable vol equal %0.32f\n',Vtot));

sys = strcat({'../lmp_mpi_'},platform,{' -in '},input,{'.in'});
system(sys{1});