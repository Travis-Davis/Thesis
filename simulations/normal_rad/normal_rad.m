% This program creates a single layer of elements for the entire simulation
% domain.  It applies a normal radius distribution to force simulation
% asymetery to force non-zero z & y force components in a pure compressive
% simulation case.

clear all
clc
system('rm atom.include *.log log.* ../data/DUMP/*.vtk ../data/DUMP/*.dump');

platform = 'mb';
input    = 'normal_rad';
%% Constants:
% Simulation Box:
x_lo = -100;
x_hi = 100;
y_lo = -100;
y_hi = 100;
z_lo = -100;
z_hi = 100;

R_o   = 1;                          % Mean Radius of particles.
sig_r = .01;                         % Standard Deviation of Particles.

% Physical:
rho_o = 900;                        % Density of particles.
rho_w = 1017;                       % Density of sea-water.
rho_a = 1;                          % Density of air.
g     = 9.81;                       % Gravitational acceleration.

Cdo   = 10;                         % Ocean Drag Coefficient.
Cdw   = 10;                          % Water Drag Coefficient.

K     = 2e5;                        % Bulk Modulus.
s00   = 1e9;                        % Bond Break criteria.
alpha = 25;                         % Bond Break criteria.

% Computed Values:
hor = 2*R_o;
c   = (18*K)/(4*hor^4);

n_x = floor((x_hi - x_lo)/(2*R_o)); % Number of Elements in x to fill sim.
n_y = floor((y_hi - y_lo)/(2*R_o)); % Number of Elements in y to fill sim.
N   = n_x*n_y;                      % Total Number of Elements.

R   = normrnd(R_o,sig_r,1,N);       % Normal Radius distribution.
V   = 4/3*pi*R.^3;                  % Volume of each element.
M   = rho_o*V;                      % Mass of each element.

%% Compute Equilibrium Position of each element:

for i=1:N
    cub = roots([rho_w*pi/3,rho_w*2*pi*R(i),-rho_w*pi*R(i)^2,-rho_o*4/3*pi*R(i)^3]);
    for j = 1:3
        if cub(j) < 0
            z_bot(i) = cub(j) - R_o;
        end
    end
end

%% Set up the Simulation:
% Simulation Constants:
fileID = fopen('constant.include','w');

% Prototype Microelastic Brittle:
fprintf(fileID,sprintf('variable c equal %0.32f\n',c));
fprintf(fileID,sprintf('variable hor equal %0.32f\n',hor));
fprintf(fileID,sprintf('variable s00 equal %0.32f\n',s00));
fprintf(fileID,sprintf('variable alpha equal %0.32f\n',alpha));

% Physical Parameters:
fprintf(fileID,sprintf('variable dens equal %0.32f\n',rho_o));
fprintf(fileID,sprintf('variable waterDensity equal %0.32f\n',rho_w));
fprintf(fileID,sprintf('variable airDensity equal %0.32f\n',rho_a));
fprintf(fileID,sprintf('variable gravity equal %0.32f\n',g));

fprintf(fileID,sprintf('variable Cdo equal %0.32f\n',Cdo));
fprintf(fileID,sprintf('variable Cdw equal %0.32f\n',Cdw));

% Create the Simulation Domain:
fileID = fopen('simbox.include','w');

fprintf(fileID,'boundary    f f f\n');
fprintf(fileID,'atom_modify map array\n');
fprintf(fileID,sprintf('variable x0 equal %0.32f\n',x_lo));
fprintf(fileID,sprintf('variable x1 equal %0.32f\n',x_hi));
fprintf(fileID,sprintf('variable y0 equal %0.32f\n',y_lo));
fprintf(fileID,sprintf('variable y1 equal %0.32f\n',y_hi));
fprintf(fileID,sprintf('variable z0 equal %0.32f\n',z_lo));
fprintf(fileID,sprintf('variable z1 equal %0.32f\n',z_hi));

fprintf(fileID,'region boxreg block ${x0} ${x1} ${y0} ${y1} ${z0} ${z1} units box\n');
fprintf(fileID,'create_box 1 boxreg');

% Particle initialization:
fileID = fopen('atom.include','w');
k = 1;
for i = 1:n_x
    for j = 1:n_y
        fprintf(fileID,sprintf('create_atoms 1 single %0.32f %0.32f %0.32f\n',...
            x_lo+(i-1)*2*R_o,y_lo+(j-1)*2*R_o,z_bot(k)));
        fprintf(fileID,sprintf('set atom %g volume %0.32f\n',...
            k,V(k)));
        fprintf(fileID,sprintf('set atom %g mass %0.32f\n',...
            k,M(k)));
        k = k+1;
    end
end
fprintf(fileID,sprintf('variable diam equal %0.32f\n',min(R)));

% Forcing initialization:

%% Run the Simulation:
sys = strcat({'../lmp_mpi_'},platform,{' -in '},input,{'.in'});
system(sys{1});

system(strcat('python ../vtk_pos_peri.py ../data/DUMP/',input,'.dump 100'));