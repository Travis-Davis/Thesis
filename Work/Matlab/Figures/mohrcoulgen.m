%% Mohr-Coulomb Envelope Generator
function [sig1plot,sig2plot,axialplot,sheerplot] = mohrcoulgen(mu_c,str_com,str_ten,res)
% Mohr-Coulomb Envelope Generator: Given a set of Axial and Sheer
% deformation limits, the function generates the set of possible failure
% stress states.  This includes both the axial/sheer stress and principal
% stress states.  The function then returns 4 vectors to the calling
% script.:
% sig1plot - principal stress state 1
% sig2plot - principal stress state 2
% axiplot  - Axial stress state
% tenplot  - Tensile stress state

%% Computation
axialplot = linspace(-str_com,str_ten,res);
sheerplot1 = -mu_c*(axialplot-str_ten);
sheerplot2 =  mu_c*(axialplot-str_ten);

sig1plot = 1/2*[(axialplot+sheerplot1),(axialplot+sheerplot2)];
sig2plot = 1/2*[(axialplot-sheerplot1),(axialplot-sheerplot2)];

axialplot = [axialplot(1),axialplot,axialplot];
sheerplot = [sheerplot2(1),sheerplot1,sheerplot2];

end