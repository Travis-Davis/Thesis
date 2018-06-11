clear all
clc

R = 1;
z = linspace(-2,2,1000);

for i = 1:length(z)
    if z(i) <= R
        if z < -2*R
            FB(i) = 4/3*pi*(R)^3;
        else
            h = z(i)-R;
            FB(i) = pi/3*((h)^3-R^2*h);
        end
    else
        FB(i) = 0;
    end
end

plot(z,FB)
