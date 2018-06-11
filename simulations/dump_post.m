clear all
clc

pos_end = length(fileread('data/DUMP/3D_trial.DUMP'));
input = fopen('data/DUMP/3D_trial.DUMP');

pos = ftell(input);
i=1;
while pos < pos_end
    Raw_Data{i} = fgetl(input);
    pos = ftell(input);
    i = i+1;
end

k = 1;
for i = 10:10:length(Raw_Data)
    handle = str2num(Raw_Data{i});

    data(k,:) = handle;

    k = k+1;
end

k = 1;
for i = 2:10:length(Raw_Data)
    handle = str2num(Raw_Data{i});

    time(k) = handle;

    k = k+1;
end

R    = .5;
Vtot = 4/3*pi*R^3;
g_r  = 1017/900*9.81;
z    = data(:,7);

for i = 1:length(z)
    if z(i) < -R
        V(i) = Vtot;
    end
    
    if z(i) > R
        V(i) = 0;
    end
    
    if z(i) >= -R && z(i) <= R
        V(i) = 1/2*Vtot+pi*(z(i)^3/3-R^2*z(i));
    end
end

%%
figure
plot(time,z)