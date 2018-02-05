clear
clc

deflim = 10^-5*3600*24*1/2*[[-1,1,-1,1];[1,-1,1,-1]];
res = 10;
P = 5E3;
ecc = 2;

sig1plot = [];
sig2plot = [];
axialplot = [];
sheerplot = [];
for i = 1:2
    [sig1plot_,sig2plot_,axialplot_,sheerplot_] = ellyieldgen(deflim(i,:),res,P,ecc);
    
    sig1plot = [sig1plot,sig1plot_];
    sig2plot = [sig2plot,sig2plot_];
    axialplot = [axialplot,axialplot_];
    sheerplot = [sheerplot,sheerplot_];
end

%% Plotting
    
figure
hold on
title({'Elliptic Yield Curve for Viscous Plastic Rheology';...
    'Principal Stress: $\sigma_1$ vs. $\sigma_2$'})
xlabel('$\sigma_1$')
ylabel('$\sigma_2$')
scatter(sig1plot,sig2plot,1,'b.')
set(gca,'visible','on')
set(gca,'XAxisLocation','origin')
set(gca,'YAxisLocation','origin')
axis([-6000,2000,-6000,2000])
plot([-P/2,-P/2],[-7000,2000],'r-')
plot([-7000,2000],[-P/2,-P/2],'r-')

figure
hold on
title({'Elliptic Yield Curve for Viscous Plastic Rheology';...
    'Axial vs. Sheer: $\sigma_I$ vs. $\sigma_{II}$'})
xlabel('$\sigma_I$')
ylabel('$\sigma_{II}$')
scatter(axialplot,sheerplot,1,'b.')
set(gca,'visible','on')
set(gca,'XAxisLocation','origin')
set(gca,'YAxisLocation','origin')
axis([-10000,1000,-3000,3000])
    