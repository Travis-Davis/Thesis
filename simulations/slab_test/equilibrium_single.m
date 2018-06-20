clear all
clc

%%
platform = 'mb';
input    = '3D_peri';

rho_w = 1017;
rho_o = 900;

R_o   = 1;

cub = roots([rho_w*pi/3,rho_w*2*pi*R_o,-rho_w*pi*R_o^2,-rho_o*4/3*pi*R_o^3]);
for i = 1:3
    if cub(i) < 0
        z_eq_ref = cub(i) - R_o;
    end
end

lim    = .5;
dist   = linspace(-lim,lim,21);
z_test = (dist+1).*z_eq_ref;
system('rm *.log log.* ../data/DUMP/*.vtk ../data/DUMP/*.dump');

%%
for i = 1:length(z_test)
    system('rm in.include.equilibrium in.trial *.log log.* ../data/DUMP/*.vtk ../data/DUMP/*.dump');
    
    fileID = fopen('in.include.equilibrium','w');
    fprintf(fileID,sprintf('create_atoms  1 single 0 0 %0.32f\n',z_test(i)));
    
    fileID = fopen('in.trial','w')
    fprintf(fileID,sprintf('variable diam equal %0.32f\n',2*R_o));
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
    
    final(2*i-1,:) = time';
    final(2*i,:) = data(:,5)' + R_o;
end

%%
figure
box on
for i = 1:length(z_test)
    hold on
    if i == floor(length(z_test)/2)+1
        plot(final(2*i-1,:)*2e-4,final(2*i,:),'LineWidth',2)
    else
        plot(final(2*i-1,:)*2e-4,final(2*i,:))
    end
end
    