clear
clc

deflim = 5000*[[-1,1,-1,1];[1,-1,1,-1]];
res = 10;

sig1plot = [];
sig2plot = [];
axialplot = [];
sheerplot = [];
for i = 1:2
    [sig1plot_,sig2plot_] = newtvisgen(deflim(i,:),res);
    
    sig1plot = [sig1plot,sig1plot_];
    sig2plot = [sig2plot,sig2plot_];
end

%% Plotting
    
figure
hold on
title({'Newtonian Fluid Rheology';...
    'Principal Stress: $\sigma_1$ vs. $\sigma_2$'})
xlabel('$\sigma_1$')
ylabel('$\sigma_2$')
scatter(sig1plot,sig2plot,1,'b.')
set(gca,'visible','on')
set(gca,'XAxisLocation','origin')
set(gca,'YAxisLocation','origin')
