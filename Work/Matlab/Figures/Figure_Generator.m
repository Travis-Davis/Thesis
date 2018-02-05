clear all
clc

%% Retrieve Data
load('NewtonianViscousRheology')
newt1 = sig1plot;
newt2 = sig2plot;

load('EllipticYieldRheology.mat')

%% Plot
figure
hold on
%title({'Elliptic Yield Curve for Viscous Plastic Rheology';...
%    'Principal Stress: $\sigma_1$ vs. $\sigma_2$'})
xlabel('$\sigma_1$')
ylabel('$\sigma_2$')
s1 = plot(sig1plot,sig2plot,'b.');
p1 = plot(newt1,newt2,'b');
plot([-P/2,-P/2],[-3000,-2000],'r-')
plot([-3000,-2000],[-P/2,-P/2],'r-')
set(gca,'visible','on')
set(gca,'XAxisLocation','origin')
set(gca,'YAxisLocation','origin')
set(gca,'PlotBoxAspectRatio',[1 1 1])
box on
legend([s1 p1],'Location','northwest','Viscous Plastic','Newtonian Viscous')
axis([-6000,3000,-6000,3000])

print('Elliptic Yield Curve.png','-dpng','-r400')