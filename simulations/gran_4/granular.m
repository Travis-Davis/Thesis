clear all
clc

%% Initialize the Run
% Declare .in file to be referenced:
input = 'gran_4';
dump  = strcat(input,'.DUMP');
atom  = 'in.include.atoms.update';
val   = 10; % number of dumped parameters

% Ice Physical Constants:
dens = 900;
mu   = 0.6;
sigt = 1e10;
sigc = -1e12;

% Model Dynamics:
% Note: Ensure these are synched with .in file
kn   = 1e6;
damp = 0.0;

% Run duration (hours):
T    = 72;

% Initial Atom Creation Region:
% Note: This MUST be in synch with the parent .in file
x0   = 0;
x1   = 500000;
y0   = x1;
y1   = x1*1.5;

% Choose Initial Floe Distribution:
%   1 = Constant Radius Lattice
%   2 = Power Law (radius) Uniform Random (spatial)
distribution = 1;

% Floe Radius Parameters:
% Note: Rmin and Rmax only needed for non-constant radius distributions
Rmean = 2000;
Rmin  = 2000;
Rmax  = 2000;

%% Run Initialization Computations
if distribution == 1
    [minmass] = unilatdist(Rmean,x0,x1,y0,y1);
end

if distribution == 2
    [minmass] = powlawdist(Rmean,Rmin,Rmax,x0,x1,y0,y1);
end

% Set the Total Run Duration:
% Note: Each run is executed for 100 timesteps
timestep = pi/sqrt(2*kn/minmass-damp^2/4.0)*0.02;
nrun = floor((T*3600)/(100*timestep));

%% Run the Simulation 
n = 1;
while n <= nrun+1
    if n == 1
        sysin = strcat('../lmp_mpi_mb -in in.i.',input);
        system(sysin)
        gran2vtk(input,n);
        clear sysin
        n = n+1;
    else
        [aind,pos,str,sca] = dump2mlab(dump,val);
        Area = pi.*(sca(2,:)).^2; 
        
        % Create a function to create a .vtk file
        if n == floor(n/5)*5
            gran2vtk(input,n);
        end
        
        % Evaluate Failure Criteria and Scale Radius
        fileID = fopen(atom,'w');
        
        fprintf(fileID,sprintf('variable minmass equal %f\n',minmass));
        fprintf(fileID,sprintf('variable diam equal %f\n',2*max(sca(2,:))));
        
        for j=1:length(aind)
            % Coulombic Failure Criteria
            foo =  str(1,j)>sigt ...
                || str(1,j)<sigc ...
                || abs(str(2,j)) > -mu*(str(1,j)-sigt);
            if foo == 1
                % This is the Deformation Relationship
                D(j) = 2*sca(2,j)*0.90;
                fprintf(fileID,sprintf('set atom %g diameter %f\n',aind(j),D(j)));
            end
        end
        
        sysin = strcat('../lmp_mpi_mb -in in.r.',input);
        system(sysin)
        n = n+1;
        clear sysin
        
        fclose(fileID);
    end
end
