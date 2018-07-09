% This program creates two multi-layer slabs of normally distributed radius
% elements.  It then allows the two slabs to reach isostatic equilibrium in
% the simulation and pushes them together to create a ridge.

clear all
clc
system('rm atom.include *.log log.* ../data/DUMP/*.vtk ../data/DUMP/*.dump');

%%
platform = 'mb';
input    = 'slab_vel';
%% Constants:
% R_o =  Mean Radius of particles; sig_r = Standard Deviation of Particles;
% v_o =  Initial Velocity of slab.
R_o   = .5;
sig_r = 0.001;
v_o   = 1.5;

% Simulation Box:
pad  = 5;

x_lo = 0;
x_hi = 30;
y_lo = 0;
y_hi = 30;
z_lo = -20;
z_hi = 20;
xc   = (x_hi+x_lo)/2;

h_1  = 1;
h_2  = 2;

x_hi_1 = xc - R_o;
x_lo_2 = xc + R_o;

% Physical:
% rho_o = Density of particles; rho_w = Density of sea-water.
% rho_a = Density of air.
rho_o = 900;
rho_w = 1017;
rho_a = 1;

% g = Gravity.
g     = 9.81;

% Cdo = Ocean drag coefficient; Cda = Air drag coefficient. 
Cdo   = 10;
Cda   = 10;

% K = Bulk Modulus; s00 = Bond Criteria; alpha = Bond Criteria.
K     = 15e3;
s00   = R_o;
alpha = 0.9;

% Computed Values:
hor   = 2*R_o;
c     = (18*K)/(4*hor^4);

% Floe 1: 
n_x_1 = round((x_hi_1-x_lo-pad)/(2*R_o));
n_y_1 = round((y_hi-y_lo-2*pad)/(2*R_o))+1;
n_z_1 = round(h_1/(2*R_o));
N_1   = n_x_1*n_y_1*n_z_1;

% Floe 2: 
n_x_2 = round((x_hi-x_lo_2)/(2*R_o))-1;
n_y_2 = round((y_hi-y_lo-2*pad)/(2*R_o))+1;
n_z_2 = round(h_2/(2*R_o));
N_2   = n_x_2*n_y_2*n_z_2;

% Element Parameters:
N = N_1+N_2;

R = normrnd(R_o,sig_r,1,N);       % Normal Radius distribution.
V = 4/3*pi*R.^3;                  % Volume of each element.
M = rho_o*V;                      % Mass of each element.

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
fprintf(fileID,sprintf('region slab1 block %0.32f %0.32f %0.32f %0.32f %0.32f %0.32f units box\n',...
                        x_lo+pad,x_hi_1,y_lo+pad,y_hi-pad,-h_1,-2*R_o));
fprintf(fileID,sprintf('region slab2 block %0.32f %0.32f %0.32f %0.32f %0.32f %0.32f units box\n',...
                        x_lo_2,x_hi,y_lo+pad,y_hi-pad,-h_2,-2*R_o));
fprintf(fileID,'create_atoms 1 region slab1\n');
fprintf(fileID,'create_atoms 1 region slab2\n');

fprintf(fileID,'group 1 region slab1\n');
fprintf(fileID,'group 2 region slab2\n');

for i = 1:N
    fprintf(fileID,sprintf('set atom %g volume %0.32f\n',...
        i,V(i)));
    fprintf(fileID,sprintf('set atom %g mass %0.32f\n',...
        i,M(i)));
end

fprintf(fileID,sprintf('variable diam equal %0.32f\n',2*min(R)));

% Slab Velocity:
fileID = fopen('velocity.include','w');
fprintf(fileID,sprintf('variable vo equal %0.32f\n',v_o));
fprintf('variable vx1 atom vx+${vo}\n');
%fprintf('variable vx2 atom vx-${vo}\n');

fprintf(fileID,'velocity 1 set ${vo} NULL NULL\n');
fprintf(fileID,'velocity 2 set -${vo} NULL NULL\n');
    
%% Run the Simulation:
sys = strcat({'../lmp_mpi_'},platform,{' -in '},input,{'.in'});
system(sys{1});

%%
system(strcat('python ../vtk_pos_peri.py ../data/DUMP/',input,'.dump 10'));