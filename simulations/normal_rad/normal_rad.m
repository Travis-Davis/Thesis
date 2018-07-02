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
x_lo = 0;
x_hi = 50;
y_lo = 0;
y_hi = 50;
z_lo = -20;
z_hi = 20;
xc   = (x_hi+x_lo)/2;

x_hi_1 = 0.95*xc;
x_lo_2 = 1.05*xc;

% R_o =  Mean Radius of particles; sig_r = Standard Deviation of Particles.
R_o   = .5;
sig_r = .001;

% Physical:
% rho_o = Density of particles; rho_w = Density of sea-water.
% rho_a = Density of air.
rho_o = 900;                        % Density of particles.
rho_w = 1017;                       % Density of sea-water.
rho_a = 1;                          % Density of air.

% g = Gravity.
g     = 9.81;

% Cdo = Ocean drag coefficient; Cda = Air drag coefficient. 
Cdo   = 10;
Cda   = 10;

% K = Bulk Modulus; s00 = Bond Criteria; alpha = Bond Criteria.
K     = 2e6;
s00   = 1e9;
alpha = .1;

% Computed Values:
hor = 2*R_o;
c   = (18*K)/(4*hor^4);

% Floe 1: 
n_x_1 = ceil((x_hi_1-x_lo)/(2*R_o));
n_y_1 = ceil((y_hi - y_lo)/(2*R_o));
N_1   = n_x_1*n_y_1;                  % Total Number of Elements in Floe 1.

% Floe 2: 
n_x_2 = ceil((x_hi-x_lo_2)/(2*R_o));
n_y_2 = ceil((y_hi - y_lo)/(2*R_o));
N_2   = n_x_2*n_y_2;                  % Total Number of Elements in Floe 2.

N = N_1+N_2;

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
fprintf(fileID,sprintf('variable Cda equal %0.32f\n',Cda));

% Create the Simulation Domain:
fileID = fopen('simbox.include','w');

fprintf(fileID,'boundary f f f\n');
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

fprintf(fileID,sprintf('lattice sc %0.32f\n',2*R_o));
fprintf(fileID,sprintf('region slab1 block %0.32f %0.32f ${y0} ${y1} %0.32f %0.32f units box\n',...
                        x_lo,x_hi_1,min(z_bot),0));
fprintf(fileID,sprintf('region slab2 block %0.32f %0.32f ${y0} ${y1} %0.32f %0.32f units box\n',...
                        x_lo_2,x_hi,min(z_bot),0));
fprintf(fileID,'create_atoms 1 region slab1\n');
fprintf(fileID,'create_atoms 1 region slab2\n');

for i = 1:N
    fprintf(fileID,sprintf('set atom %g volume %0.32f\n',...
        i,V(i)));
    fprintf(fileID,sprintf('set atom %g mass %0.32f\n',...
        i,M(i)));
end

fprintf(fileID,sprintf('variable diam equal %0.32f\n',2*min(R)));

%% Run the Simulation:
sys = strcat({'../lmp_mpi_'},platform,{' -in '},input,{'.in'});
system(sys{1});

system(strcat('python ../vtk_pos_peri.py ../data/DUMP/',input,'.dump 100'));