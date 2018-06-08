clear
clc

%% Constants
deflim = 10^-5*3600*24*1/2*[[-1,1,-1,1];[1,-1,1,-1]];
res = 20;
P = 5E3;
ecc = 2;

mu_c = tand(15);
str_com = 2*P;
str_ten = P/10;

%% Viscous-Plastic Rheology Gen.
sig1plot_VP = [];
sig2plot_VP = [];
axialplot_VP = [];
sheerplot_VP = [];
for i = 1:2
    [sig1plot_,sig2plot_,axialplot_,sheerplot_] = ellyieldgen(deflim(i,:),res,P,ecc);
    
    sig1plot_VP  = [sig1plot_VP,sig1plot_];
    sig2plot_VP  = [sig2plot_VP,sig2plot_];
    axialplot_VP = [axialplot_VP,axialplot_];
    sheerplot_VP = [sheerplot_VP,sheerplot_];
end

pos_sheer_vp = sheerplot_VP(sheerplot_VP >= 0);
pos_axial_vp = axialplot_VP(sheerplot_VP >= 0);
neg_sheer_vp = sheerplot_VP(sheerplot_VP < 0);
neg_axial_vp = axialplot_VP(sheerplot_VP < 0);

[pos_axial_vp,pos] = sort(pos_axial_vp,'ascend');
[neg_axial_vp,neg] = sort(neg_axial_vp,'descend');

axialplot_VP = [pos_axial_vp,neg_axial_vp];
sheerplot_VP = [pos_sheer_vp(pos),neg_sheer_vp(neg)];

%% Mohr-Coulomb Rheology Gen.
[sig1plot_mc,sig2plot_mc,axialplot_mc,sheerplot_mc] = mohrcoulgen(mu_c,str_com,str_ten,res);

%% Plotting
pad = 500;
figpad = .005;

xmin = min([min(axialplot_VP),min(axialplot_mc)])-pad;
xmax = max([max(axialplot_VP),max(axialplot_mc)])+pad;
ymin = min([min(sheerplot_VP),min(sheerplot_mc)])-pad;
ymax = max([max(sheerplot_VP),max(sheerplot_mc)])+pad;

figure(1)
hold on

s1 = plot(axialplot_VP,sheerplot_VP,'Color','b','LineWidth',2.5);
s2 = plot(axialplot_mc,sheerplot_mc,'Color','k','LineWidth',2.5);

xlabel('$\mathbf{\sigma_{I}}$','FontSize',14)
ylabel('$\mathbf{\sigma_{II}}$','FontSize',14)
axis([xmin,xmax,ymin,ymax])

lgd = legend([s1 s2],'Location','northwest','Viscous Plastic','Mohr-Coulomb');
lgd.Interpreter = 'latex';

fig = gcf; ax = gca;
figureformat(fig,ax,figpad)

print(fig,'Elliptic_Mohr_Coulomb','-depsc')