clear
clc

%% Test Case
lat = [80,60,-60];
lon = [.1,45, 45];
SGN = [ 1, 1, -1];

x = [];
y = [];
latr = [];
lonr = [];

for i = 1:length(lat)
    [x(i),y(i)]       = ncgeodetictoxy(lat(i),lon(i),SGN(i));
    [latr(i),lonr(i)] = ncxytogeodetic(  x(i),  y(i),SGN(i));
end

testlat = lat./latr;
testlon = lon./lonr;

disp(testlat)
disp(testlon)

