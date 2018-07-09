clear all
clc

x = [0,1,2];
y = [0,1,2];
z = [0,1,2];

n = 1;
for i = 1:length(x)
    for j = 1:length(y)
        for k = 1:length(z)
            vol(n,:) = [x(i),y(j),z(k),0];
            n = n+1;
        end
    end
end

% Not needed if bond list is available:
hor = 1.5;
neigh = {};

for i = 1:length(vol)
    place = i-1;
    for j = 1:length(vol)
        d=(vol(i,1)-vol(j,1))^2+(vol(i,2)-vol(j,2))^2+(vol(i,3)-vol(j,3))^2;
        if d < hor && d ~= 0
            vol(i,4) = vol(i,4)+1;
            place(end+1) = j-1;
        end
        neigh{i} = place;
    end
end

DT = delaunayTriangulation(vol(:,1),vol(:,2),vol(:,3));
figure
tetramesh(DT,'FaceAlpha',0.3);

[K,v] = convexHull(DT);
figure
trisurf(K,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3))
