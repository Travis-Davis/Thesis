clear all
clc

%%
platform = 'mb';
input    = '3D_peri';

tol = 1e-5;

rho_w = 1017;
rho_o = 900;

R_o   = 1;

cub = roots([1,-3*R_o,0,4*(rho_w-rho_o)/rho_w*R_o^3]);
for i = 1:3
    if cub(i) == abs(cub(i))
        if cub(i) < R_o && cub(i) > 0
            z_eq_ref = cub(i)-R_o;
        end
    end
end

lim = 2;
dist = linspace(-lim,lim,11);
z_test = (dist+1).*z_eq_ref;
system('rm *.log log.* ../data/DUMP/*.vtk ../data/DUMP/*.dump');

%%
variance = 1;
z_eq = z_eq_ref;
z_eq_trend = z_eq_ref;

while variance > tol
    system('rm in.include.equilibrium in.trial *.log log.* ../data/DUMP/*.vtk ../data/DUMP/*.dump');
    
    fileID = fopen('in.include.equilibrium','w');
    fprintf(fileID,sprintf('create_atoms  1 single 0 0 %0.32f\n',z_eq));
    
    fileID = fopen('in.trial','w')
    fprintf(fileID,sprintf('variable diam equal %f\n',2*R_o));
    fprintf(fileID,sprintf('variable rad equal ${diam}*0.5\n'));
    
    system(strcat('../lmp_mpi_',platform,' -in in.',input));
    
    pos_end = length(fileread(strcat('../data/DUMP/',input,'.dump')));
    dump = fopen(strcat('../data/DUMP/',input,'.dump'));

    pos = ftell(dump);
    
    j=1;
    while pos < pos_end
        Raw_Data{j} = fgetl(dump);
        pos = ftell(dump);
        j = j+1;
    end
    
    j = 1;
    for k = 10:10:length(Raw_Data)
        handle = str2num(Raw_Data{k});

        data(j,:) = handle;

        j = j+1;
    end
    
    j = 1;
    for k = 2:10:length(Raw_Data)
        handle = str2num(Raw_Data{k});

        time(j) = handle;

        j = j+1;
    end
    
    final(1,:) = time';
    final(2,:) = data(:,5)';
    variance = var(data(:,5)');
    
    del  = final(2,2)-final(2,1);
    z_eq = z_eq+del;
    z_eq_trend(end+1) = z_eq;
end

%%
figure
box on
plot(final(1,:),final(2,:),'k')