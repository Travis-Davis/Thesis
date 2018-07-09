clear all
clc
%%
R_o   = .5;
rho_o = 900;
sig_r = 0.001;
h1   = 2;
h2   = 3;

pad  = 5;

x_lo = 0;
x_hi = 30;
y_lo = 0;
y_hi = 30;
z_lo = -20;
z_hi = 20;
xc   = (x_hi+x_lo)/2;


x1_lo = x_lo+pad;
x1_hi = xc-2*R_o;
y1_lo = y_lo+pad;
y1_hi = y_hi-pad;
z1_lo = -h1-2*R_o;
z1_hi = -2*R_o;

n_x_1 = round((x1_hi-x1_lo)/(2*R_o));
n_y_1 = round((y1_hi-y1_lo)/(2*R_o))+1;
n_z_1 = round(h1/(2*R_o));
N_1   = n_x_1*n_y_1*n_z_1;

x2_lo = xc+2*R_o;
x2_hi = x_hi-pad;
y2_lo = y_lo+pad;
y2_hi = y_hi-pad;
z2_lo = -h2-2*R_o;
z2_hi = -2*R_o;

n_x_2 = round((x2_hi-x2_lo)/(2*R_o))-1;
n_y_2 = round((y2_hi-y2_lo)/(2*R_o))+1;
n_z_2 = round(h2/(2*R_o));
N_2   = n_x_2*n_y_2*n_z_2;

x1 = linspace(x1_lo,x1_hi+2*R_o,n_x_1);
y1 = linspace(y1_lo,y1_hi+2*R_o,n_y_1);
z1 = linspace(z1_lo,z1_hi+2*R_o,n_z_1);

x2 = linspace(x2_lo,x2_hi+2*R_o,n_x_2);
y2 = linspace(y2_lo,y2_hi+2*R_o,n_y_2);
z2 = linspace(z2_lo,z2_hi+2*R_o,n_z_2);

ind1 = 0;
for i = 1:n_z_1
    for j = 1:n_y_1
        for k = 1:n_x_1
            ind1 = ind1+1;
            p1(ind1,1) = ind1-1;
            p1(ind1,2) = x1(k);
            p1(ind1,3) = y1(j);
            p1(ind1,4) = z1(i);
        end
    end
end

ind2 = 0;
for i = 1:n_z_2
    for j = 1:n_y_2
        for k = 1:n_x_2
            ind2 = ind2+1;
            p2(ind2,1) = ind2-1+N_1;
            p2(ind2,2) = x2(k);
            p2(ind2,3) = y2(j);
            p2(ind2,4) = z2(i);
        end
    end
end

DT1 = delaunayTriangulation(p1(:,2),p1(:,3),p1(:,4));
DT2 = delaunayTriangulation(p2(:,2),p2(:,3),p2(:,4));

den = normrnd(rho_o,sig_r,1,N_1+N_2);
x1 = DT1.Points(:,1);
x2 = DT2.Points(:,1);
y1 = DT1.Points(:,2);
y2 = DT2.Points(:,2);
z1 = DT1.Points(:,3);
z2 = DT2.Points(:,3);
C1 = DT1.ConnectivityList;
C2 = DT2.ConnectivityList;

cellnum = length(C1(:,1))+length(C2(:,1));

%%
vtk = fopen('test.vtk','w');

fprintf(vtk,'# vtk DataFile Version 3.0\n');
fprintf(vtk,'DUMP output\n');
fprintf(vtk,'ASCII\n');
fprintf(vtk,'DATASET UNSTRUCTURED_GRID\n');
fprintf(vtk,sprintf('POINTS %g float\n',N_1+N_2));

for i = 1:N_1
    fprintf(vtk,sprintf('%0.8f %0.8f %0.8f\n',x1(i),y1(i),z1(i)));
end

for i = 1:N_2
    fprintf(vtk,sprintf('%0.8f %0.8f %0.8f\n',x2(i),y2(i),z2(i)));
end

fprintf(vtk,sprintf('CELLS %g %g\n',cellnum,5*cellnum));
nc = 1;
for i = 1:length(C1)
    fprintf(vtk,sprintf('4 %g %g %g %g\n',C1(i,:)-1));
    nc = nc+1;
end
for i = 1:length(C2)
    fprintf(vtk,sprintf('4 %g %g %g %g\n',C2(i,:)+N_1-1));
    nc = nc+1;
end

fprintf(vtk,sprintf('CELL_TYPES %g\n',nc-1));
for i = 1:nc-1
    fprintf(vtk,'10\n');
end

fprintf(vtk,sprintf('POINT_DATA %g\n',N_1+N_2));
fprintf(vtk,'SCALARS density float 1\n');
fprintf(vtk,'LOOKUP_TABLE default\n');
for i=1:length(den)
    fprintf(vtk,sprintf('%0.8f\n',den(i)));
end