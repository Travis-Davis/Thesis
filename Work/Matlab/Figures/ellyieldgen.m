%% Elliptic Yield Curve Generator
function [sig1plot,sig2plot,axialplot,sheerplot] = ellyieldgen(deflim,res,P,ecc)
% Elliptic Yield Curve Generator: Given a set of Axial and Sheer
% deformation limits, the function generates the set of possible stress
% states.  This includes both the axial/sheer stress and principal stress
% states.  The function then returns 4 vectors to the calling script.

%% Constants
axialmin = deflim(1);
axialmax = deflim(2);
sheermin = deflim(3);
sheermax = deflim(4);
n = res;

dux = linspace(axialmin,axialmax,n);
duy = linspace(sheermin,sheermax,n);
dvx = linspace(sheermin,sheermax,n);
dvy = linspace(axialmin,axialmax,n);

%% Pre-Allocation
% Strain Space
str = {};
strt = {};

strt{n,n} = [];
str{n,n} = strt;

% Stress Space
sig = {};
sigt = {};

sigt{n,n} = [];
sig{n,n} = sigt;

% Principal Stress Space
sigp = {};
sigpt = {};

sigpt{n,2*n} = [];
sigp{n,2*n} = sigpt;

% Non-Linear Shear Viscosity
eta = {};

eta{n,n} = [];

% Non-Linear Bulk Viscosity
zeta = {};

zeta{n,n} = [];

%% Define Strain Rate, Bulk Strength, Sheer Strength, and Stress spaces
for i = 1:n
    for j = 1:n
        for k = 1:n
            for l = 1:n
                e11 = dux(i);
                e22 = dvy(l);
                e12 = 1/2*(dvx(k)+duy(j));
                e21 = 1/2*(dvx(k)+duy(j));
                
                del_ = sqrt((e11^2+e22^2)*(1+1/ecc^2)+4/ecc^2*e12^2+...
                    2*e11*e22*(1-1/ecc^2));
                zeta_ = P/(2*del_);
                eta_ = P/(2*del_*ecc^2);
                
                sig11 = 2*eta_*e11+(zeta_-eta_)*(e11+e22)-P/2;
                sig12 = 2*eta_*e12;
                sig21 = 2*eta_*e21;
                sig22 = 2*eta_*e22+(zeta_-eta_)*(e11+e22)-P/2;
                det = sig11*sig22 - sig12*sig21;
                
                sigp1a = 1/2*(sig11+sig22+sqrt((sig11+sig22)^2-4*det));
                sigp2a = 1/2*(sig11+sig22-sqrt((sig11+sig22)^2-4*det));
                sigp1b = 1/2*(sig11+sig22-sqrt((sig11+sig22)^2-4*det));
                sigp2b = 1/2*(sig11+sig22+sqrt((sig11+sig22)^2-4*det));
                
                str{i,j}{k,l}(1,1) = e11;
                str{i,j}{k,l}(1,2) = e12;
                str{i,j}{k,l}(2,1) = e21;
                str{i,j}{k,l}(2,2) = e22;
                
                del{i,j}(k,l) = del_;
                zeta{i,j}(k,l) = zeta_;
                eta{i,j}(k,l) = eta_;
            
                sig{i,j}{k,l}(1,1) = sig11;
                sig{i,j}{k,l}(1,2) = sig12;
                sig{i,j}{k,l}(2,1) = sig21;
                sig{i,j}{k,l}(2,2) = sig22;
            
                sigp{i,j}{k,l}(1,1) = sigp1a;
                sigp{i,j+n}{k,l}(1,1) = sigp1b;
                sigp{i,j}{k,l}(1,2) = 0;
                sigp{i,j+n}{k,l}(1,2) = 0;
                sigp{i,j}{k,l}(2,1) = 0;
                sigp{i,j+n}{k,l}(2,1) = 0;
                sigp{i,j}{k,l}(2,2) = sigp2a;
                sigp{i,j+n}{k,l}(2,2) = sigp2b;
                      
                axial(k,l) = (sig11+sig22);
                sheer(k,l) = sqrt((sig11+sig22)^2-4*(sig11*sig22-sig12*sig21));
                sheer(k,l+n) = -sqrt((sig11+sig22)^2-4*(sig11*sig22-sig12*sig21));
            end
        end
    end
end

%% Output Transformations
m = 1;
for i = 1:n
    for j = 1:2*n
        for k = 1:n
            for l = 1:n
                sig1plot(m) = sigp{i,j}{k,l}(1,1);
                sig2plot(m) = sigp{i,j}{k,l}(2,2);
                
                m = m+1;
            end
        end
    end
end

axialplot = [];
sheerplot = [];
for i = 1:n
    axialplot = [axialplot,axial(i,:),axial(i,:)];
    sheerplot = [sheerplot,sheer(i,:)];
end
end