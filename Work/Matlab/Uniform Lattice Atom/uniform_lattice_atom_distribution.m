clear all
clc

%% Lattice Spatial Distribution//Uniform Radius Distribution
R   = 2000;
rho = 900;
h   = 1;

mass = rho*h*pi*R^2;
xDIM = 500000;
yDIM = 1.5*xDIM;

% Atom Region Specifications
x0 = 100000;
x1 = xDIM-x0;
y0 = xDIM;
y1 = yDIM;

nx = floor((x1-x0)/(2*R));
ny = floor((y1-y0)/(2*R));

N = nx*ny;

%% x,y Coordinate Generation
n=1;
for i=1:nx
    for j=1:ny
        x(n) = x0+(2*i-1)*R;
        y(n) = y0+(2*j-1)*R;
        
        n=n+1;
    end
end


%% .in Script Generation
% Element Generation
fileID = fopen('in.include.atoms.grid','w');

fprintf(fileID,'###################\n');
fprintf(fileID,'# Atom Parameters #\n');
fprintf(fileID,'###################\n\n');

fprintf(fileID,sprintf('variable minmass equal %f\n',mass));
fprintf(fileID,sprintf('variable diam equal %f\n',2*R));

for i =  1:N
    fprintf(fileID,sprintf('create_atoms 1 single %f %f %f\n',x(i),y(i),0));
end

for i =  1:N
    fprintf(fileID,sprintf('set   atom %g diameter %f\n',i,2*R));
    fprintf(fileID,sprintf('set   atom %g mass %f\n',i,mass));
end