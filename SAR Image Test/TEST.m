clear all
clc

input = 'S1B_ESA_2018_01_22_11_09_12_0569934552_030.25W_84.16N_HH_C5_GFS05CDF_wind.tif';

meta = geotiffinfo(input);
I = imread(input);
latlim = meta.SpatialRef.LatitudeLimits;
lonlim = meta.SpatialRef.LongitudeLimits;

I = rgb2gray(I);
edge = edge(I,'Prewitt');
latedge = linspace(latlim(1),latlim(2),length(edge(:,1)));
lonedge = linspace(lonlim(1),lonlim(2),length(edge(1,:)));

sample = 1;
[Z,refvec]  = etopo(getenv('ETOPOFILE'),sample,latlim,lonlim);
Z = flipud(Z);

latz = linspace(latlim(1),latlim(2),length(Z(:,1)));
lonz = linspace(lonlim(1),lonlim(2),length(Z(1,:)));

for i = 1:length(edge(:,1))
    for j = 1:length(edge(1,:))
        latvec = abs(latz - latedge(i));
        lonvec = abs(lonz - lonedge(j));
        
        [zlat,latind] = min(latvec);
        [zlon,lonind] = min(lonvec);
        
        land(i,j) = Z(latind,lonind);
    end
end

land = land > 0;
output = land | edge;
imshow(output)