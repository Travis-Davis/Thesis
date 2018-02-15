clear all
clc

%% Power-Law Ice Floe Generator
%
%
%% Parameters
N    = 5000;
mean = 1250;
rmin = 500;
rmax = 20000;
rho  = 900;
h    = 1;

% Atom Region Specifications
x0 = 0;
x1 = 500000;
y0 = x1;
y1 = x1*1.5;

%% Generate Elements
i = 1;
while i <= N
    radt = random('Exponential',mean);
    if radt > rmin && radt < rmax
        rad(i) = radt;
        i = i+1;
    end
end

rad = sort(rad,'descend');
mass = rho*h*pi*(rad).^2;

diam = max(rad);
minmass = min(mass);

area = sum(pi*(rad).^2);
domain = (y1-y0)*(x1-x0);

if 1.25*area > domain
    fprintf('Overpack: Reduce N or r,mean\n')
    return
end

%% Generate x,y coordinates

res = min(rad)/2;
x_n_mask = floor((x1-x0)/res);
y_n_mask = floor((y1-y0)/res);

x = linspace(x0,x1,x_n_mask);
y = linspace(y0,y1,y_n_mask);
mask = zeros(y_n_mask,x_n_mask);

k = 1;
while k <= N
    temp_mask = zeros(y_n_mask,x_n_mask);
    t_dim = ceil(rad(k)/res);
    
    l = 1;
    while l == 1
        x_draw = random('Uniform',x0,x1);
        y_draw = random('Uniform',y0,y1);
        
        x_diff = abs(x-x_draw);
        y_diff = abs(y-y_draw);
        x_loc = find(x_diff == min(x_diff));
        y_loc = find(y_diff == min(y_diff));
        
        if x_loc-t_dim > 0 && x_loc+t_dim < length(x)
            if y_loc-t_dim > 0 && y_loc+t_dim < length(y)
                l = l+1;
            end
        end
    end

    for i=1:2*t_dim
        for j=1:2*t_dim
            pos = sqrt((i-t_dim)^2+(j-t_dim)^2);
            if pos <= t_dim
                temp_mask(y_loc-t_dim+i,x_loc-t_dim+j) = 1;
            else
                temp_mask(y_loc-t_dim+i,x_loc-t_dim+j) = 0;
            end
        end
    end

    bool_mask = temp_mask + mask;
    
    if any(bool_mask(:) >= 2)
    else
        mask = bool_mask;
                
        x_ele(k) = x(x_loc);
        y_ele(k) = y(y_loc);
        
        k = k+1;
    end

end

%% .in Script Generation
% Element Generation
fileID = fopen('in.include.atoms.power','w');

fprintf(fileID,'###################\n');
fprintf(fileID,'# Atom Parameters #\n');
fprintf(fileID,'###################\n\n');

fprintf(fileID,sprintf('variable minmass equal %f\n',minmass));
fprintf(fileID,sprintf('variable diam equal %f\n',diam));

for i =  1:N
    fprintf(fileID,sprintf('create_atoms 1 single %f %f %f\n',x_ele(i),y_ele(i),0));
end

for i =  1:N
    fprintf(fileID,sprintf('set   atom %g diameter %f\n',i,2*rad(i)));
    fprintf(fileID,sprintf('set   atom %g mass %f\n',i,mass(i)));
end