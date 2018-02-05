%% SAR Image Edge Detector (Prewitt Edge Finding Algorithm)
function [output] = sarprewitt(input)
% Raw NRCS .tiff file ingest and processing
meta = geotiffinfo(input);
I = imread(input);
latlim = [meta.BoundingBox(1,2),meta.BoundingBox(2,2)];
lonlim = [meta.BoundingBox(1,1),meta.BoundingBox(2,1)];

I = rgb2gray(I);
edge = edge(I,'Prewitt');
latedge = linspace(latlim(1),latlim(2),length(edge(:,1)));
lonedge = linspace(lonlim(1),lonlim(2),length(edge(1,:)));

sample = 1;
[Z,refvec]  = etopo(getenv('ETOPOFILE'),sample,latlim,lonlim);
latz = linspace(latlim(1),latlim(2),length(Z(1,:)));
lonz = linspace(latlim(1),latlim(2),length(Z(:,1)));

for i = 1:length(edge(:,1))
    for j = 1:length(edge(1,:))
        latvec = abs(latz - latedge(i));
        lonvec = abs(lonz - lonedge(j));
        
        zlat = latz(latvec == min(latvec));
        zlon = lonz(latvec == min(lonvec));
        
        land(i,j) = Z(zlat,zlon);
        
        clear latvec lonvec zlat zlon
    end
end

Z = Z' > 0;

RGB = imread(input);
I = rgb2gray(RGB);
BW = edge(I,'Prewitt');

output = Z | BW;


