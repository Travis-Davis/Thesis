function [coast] = coastgen[latmin,latmax,lonmin,lonmax]
landm = shaperead('landareas.shp','UseGeoCoords',true);
lon_all = [landm.Lon];
lat_all = [landm.Lat];

lon_region = lon_all(lon_all <= lon_max && lon_all >= lon_min);
lat_region = lat_all(lat_all <= lat_max && lat_all >= lat_max);

latlim = [latmin,latmax];
lonlim = [lonmin,lonmax];
sample = 1;
[Z,refc] = etopo(getenv('ETOPOFILE'),sample,latlim,lonlim);

if Z >= 0
    

%% Constants
num_particles = 1; %10000 particles with mean diameter of 50km
pp_seg = 5;
dia_mean = 50000;
density = 900;
lat_max = 65;
lat_circ = 66.5622;
[rx_circ,place] = ncgeodetictoxy(lat_circ,90,1);
[place,ry_circ] = ncgeodetictoxy(lat_circ,180,1);
clear place

r_circ = sqrt(rx_circ^2+ry_circ^2);

%% Coast line Element Generation
landm = shaperead('landareas.shp','UseGeoCoords',true);
lon_all = [landm.Lon];
lat_all = [landm.Lat];

lon_region = lon_all(lat_all > lat_max);
lat_region = lat_all(lat_all > lat_max);

x_arctic = [];
y_arctic = [];
del      = [];

for i = 1:length(lat_region)
    [x_arctic(i),y_arctic(i)] = ncgeodetictoxy(lat_region(i),lon_region(i),1);
end

for i = 1:length(x_arctic)-1
    del(i) = sqrt((x_arctic(i+1)-x_arctic(i))^2 + ...
                  (y_arctic(i+1)-y_arctic(i))^2);
end

vert_dia = min(del(del ~= 0));

x_arc = [];
y_arc = [];
dia_arc = [];

for i = 1:length(del)
     x_arc(end+1) = x_arctic(i);
     y_arc(end+1) = y_arctic(i);
     dia_arc(end+1) = vert_dia;
     
     if del(i) < 90 %% Value derived from Histogram of line segment lengths
        l_seg = del(i) - vert_dia;
        
        if l_seg > 0
            theta = atand((y_arctic(i+1) - y_arctic(i))/(x_arctic(i+1) - x_arctic(i)));
            
            if x_arctic(i+1) - x_arctic(i) < 0
                theta = theta + 180;
            end

            seg_dia = l_seg/pp_seg;
            
            for j = 1:pp_seg
                del_seg = 1/2*(vert_dia+seg_dia)+(j-1)*seg_dia;

                x_arc(end+1) = x_arctic(i) + cosd(theta)*del_seg;
                y_arc(end+1) = y_arctic(i) + sind(theta)*del_seg;
                dia_arc(end+1) = seg_dia;
            end
        end
     end
end

%% Z-Filter Generation
latlim = [60,90];
lonlim = [-180,180];
sample = 1;
[Z,refc] = etopo(getenv('ETOPOFILE'),sample,latlim,lonlim);

lat_vec = linspace(latlim(1),latlim(2),length(Z(:,1)));
lon_vec = linspace(lonlim(1),lonlim(2),length(Z(1,:)));

z_bool = [];

for i = 1:length(lat_vec)
    for j = 1:length(lon_vec)
        if Z(i,j) > 0
            z_bool(i,j) = 0;
        else
            z_bool(i,j) = 1;
        end
    end
end

%% Ice Element Generation
n = (1:num_particles)';
diam  = [];
x_ele = [];
y_ele = [];
z_ele = [];

lat_draw = lat_circ-1;

k = 1;

while k ~= num_particles+1;
    diam(k) = random('uniform',dia_mean*1/2,dia_mean*3/2);
    
    while lat_draw < lat_circ
        x_draw = random ('uniform',-rx_circ,rx_circ);
        y_draw = random ('uniform',-ry_circ,ry_circ);
        
        [lat_draw,lon_draw] = ncxytogeodetic(x_draw,y_draw,1);
    end
    
    [lat_ele,lon_ele] = ncxytogeodetic(x_draw,y_draw,1);
    
    lat_bool = abs(lat_vec - lat_ele);
    lon_bool = abs(lon_vec - lon_ele);
    
    [place,z_lat_ind] = min(lat_bool);
    [place,z_lon_ind] = min(lon_bool);
    
    if z_bool(z_lat_ind,z_lon_ind) == 1
        z_ele(k) = 0;
        x_ele(k) = 1000*x_draw;
        y_ele(k) = 1000*y_draw;
        k = k+1;
    end
    
    lat_draw = lat_circ - 1;
end

mass = pi*(diam/2).^2*density;
zlim = max(diam)/2;
dia = max(diam);
minmass = min(mass);

%% Transform to LAMMPS Coordinates (Lower Left = 0,0)
x_arctic = x_arctic * 1000;
y_arctic = y_arctic * 1000;
x_arc = x_arc * 1000;
y_arc = y_arc * 1000;
dia_arc = dia_arc * 1000;

x_ele = x_ele;
y_ele = y_ele;

lx = 2*(max([abs(x_ele),abs(x_arctic)])+100000);
ly = 2*(max([abs(y_ele),abs(y_arctic)])+100000);

x_arctic = x_arctic + lx;
y_arctic = y_arctic + ly;
x_arc = x_arc + lx;
y_arc = y_arc + ly;
x_ele = x_ele + lx;
y_ele = y_ele + ly;
%% Validation
for i = 1:3600
    [x_circ(i),y_circ(i)] = ncgeodetictoxy(66.5622,i/10,1);
end

x_circ = x_circ*1000 + lx;
y_circ = y_circ*1000 + ly;

figure
plot(x_circ,y_circ,'r')
hold on
for i = 1:length(del)
     if del(i) < 90
        plot([x_arctic(i),x_arctic(i+1)],...
             [y_arctic(i),y_arctic(i+1)],'b')
    end
end

scatter(x_ele,y_ele)
axis([lx/2,3/2*lx,ly/2,3/2*ly])

figure
hist(log(del))

figure
plot(x_circ,y_circ,'r')
hold on
sz = 12;
scatter(x_arc,y_arc,sz,'.r')
axis([0,2*lx,0,2*ly])
title('Simulation Boundaries')

%% .in Script Generation
% Element Generation
fileID = fopen('atoms.in','w');

fprintf(fileID,'###################\n');
fprintf(fileID,'# Atom Parameters #\n');
fprintf(fileID,'###################\n\n');

fprintf(fileID,sprintf('variable minmass equal %f\n',minmass));
fprintf(fileID,sprintf('variable dia equal %f\n',dia));

for i =  1:length(x_ele)
    fprintf(fileID,sprintf('create_atoms 1 single %f %f %f\n',x_ele(i),y_ele(i),0));
end

for i =  1:num_particles
    fprintf(fileID,sprintf('set   atom %g diameter %f\n',n(i)+length(x_arc),diam(i)));
    fprintf(fileID,sprintf('set   atom %g mass %f\n',n(i)+length(x_arc),mass(i)));
end

fprintf(fileID,'\ncompute rad all property/atom radius\n');
fprintf(fileID,'variable area atom PI*c_rad^2');

% Boundary Generation
fileID = fopen('boundary.in','w');

fprintf(fileID,'#######################\n');
fprintf(fileID,'# Boundary Parameters #\n');
fprintf(fileID,'#######################\n\n');

fprintf(fileID,'boundary f f p\n\n');

fprintf(fileID,'variable x0 equal 0\n');
fprintf(fileID,sprintf('variable x1 equal %f\n',2*lx));
fprintf(fileID,'variable y0 equal 0\n');
fprintf(fileID,sprintf('variable y1 equal %f\n',2*ly));
fprintf(fileID,sprintf('variable zlim equal %f\n\n',zlim));

fprintf(fileID,'region simbox block ${x0} ${x1} ${y0} ${y1} -${zlim} ${zlim}\n');

fprintf(fileID,'create_box 1 simbox\n\n');

for i =  1:length(x_arc)
    fprintf(fileID,sprintf('create_atoms 1 single %f %f %f\n',x_arc(i),y_arc(i),0));
end

for i =  1:length(x_arc)
    fprintf(fileID,sprintf('set   atom %g diameter %f\n',i,dia_arc(i)));
end
fprintf(fileID,'group coast region simbox\n');