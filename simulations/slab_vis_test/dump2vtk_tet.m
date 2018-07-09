function dump2vtk_tet(name,N_1)
% Read in .dump file:
str = strcat('../data/DUMP/',name,'.dump');

pos_end = length(fileread(str));
input = fopen(str);

pos = ftell(input);
i=1;
while pos < pos_end
    Raw_Data{i} = fgetl(input);
    pos = ftell(input);
    i = i+1;
end

% Scan The file for 'ITEM: TIMESTEP':
% This is done to find the time step breaks.
b_line = [];
for i = 1:length(Raw_Data)
    if length(Raw_Data{i})==14 && Raw_Data{i}(7)=='T'
        b_line(end+1) = i;
    end
end

% Create the initial 3D Triangulation:
for i = 10:b_line(2)-1
    pts(i-9,:) = str2num(Raw_Data{i});
end
pts = sortrows(pts);

p1  = pts(1:N_1,:);
p2  = pts(N_1+1:end,:);

DT1 = delaunayTriangulation(p1(:,2),p1(:,3),p1(:,4));
DT2 = delaunayTriangulation(p2(:,2),p2(:,3),p2(:,4));

x1 = DT1.Points(:,1);
x2 = DT2.Points(:,1);
y1 = DT1.Points(:,2);
y2 = DT2.Points(:,2);
z1 = DT1.Points(:,3);
z2 = DT2.Points(:,3);
C1 = DT1.ConnectivityList;
C2 = DT2.ConnectivityList;

cellnum = length(C1(:,1))+length(C2(:,1));

% Write .vtk files:
for i = 1:length(b_line)
    n = 1;
    if i == length(b_line)
        for j = b_line(i)+9:length(Raw_Data)
            pts(n,:) = str2num(Raw_Data{j});
            n=n+1;
        end
    else
        for j = b_line(i)+9:b_line(i+1)-1
            pts(n,:) = str2num(Raw_Data{j});
            n=n+1;
        end
    end
    pts = sortrows(pts);
    
    vtk = fopen(strcat('../data/DUMP/',name,num2str(i-1),'.vtk'),'w');
    
    fprintf(vtk,'# vtk DataFile Version 3.0\n');
    tst = strcat({'Time step '},Raw_Data{b_line(i)+1},{'\n'});
    fprintf(vtk,tst{1});
    fprintf(vtk,'ASCII\n');
    fprintf(vtk,'DATASET UNSTRUCTURED_GRID\n');
    fprintf(vtk,sprintf('POINTS %g float\n',length(pts(:,1))));

    for j = 1:length(pts(:,1))
        fprintf(vtk,sprintf('%0.8f %0.8f %0.8f\n',pts(j,[2,3,4])));
    end

    fprintf(vtk,sprintf('CELLS %g %g\n',cellnum,5*cellnum));

    for j = 1:length(C1)
        fprintf(vtk,sprintf('4 %g %g %g %g\n',C1(j,:)-1));
    end
    for j = 1:length(C2)
        fprintf(vtk,sprintf('4 %g %g %g %g\n',C2(j,:)+N_1-1));
    end

    fprintf(vtk,sprintf('CELL_TYPES %g\n',cellnum));
    for j = 1:cellnum
        fprintf(vtk,'10\n');
    end

    fprintf(vtk,sprintf('POINT_DATA %g\n',length(pts(:,1))));
    fprintf(vtk,'SCALARS density float 1\n');
    fprintf(vtk,'LOOKUP_TABLE default\n');
    fprintf(vtk,sprintf('%0.8f\n',pts(:,5)));

    fprintf(vtk,'SCALARS damage float 1\n');
    fprintf(vtk,'LOOKUP_TABLE default\n');
    fprintf(vtk,sprintf('%0.8f\n',pts(:,6)));
    
    disp(i-1)
    fclose(vtk);
end

end