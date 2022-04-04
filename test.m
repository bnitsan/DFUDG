clear
addpath(genpath('src'))

%%
HaloParamBar.type = 1; % 1-Sersic
HaloParamBar.nSersic = 0.61;
HaloParamBar.MSersic = 1.2e3;
HaloParamBar.rSersic = 1.9;
HaloParamBar.rStall = (HaloParamBar.rSersic)*0.3;
HaloParameters = HaloParamBar;

[rS,densityS,massS,sigmaS,phiS]=getSersic(HaloParameters.nSersic,HaloParameters.MSersic,HaloParameters.rSersic);

rc = 1:0.1:50;
rc = rc';
sigma=interp1(rS,sigmaS,rc);

%%
clear
[r,Density,Mass,Sigma,Phi] = getSersic(0.61,1e2,1.9);
[r2,Density2,Mass2,Sigma2,Phi2] = getSersic(2.0,1e2,1.9);

figure
plot(r,Density)
hold on
plot(r2,Density2)

grid on;
set(gca,'FontSize',15);
set(gca,'XScale','log')
set(gca,'YScale','log')

figure
plot(r,Mass/Mass(end))
hold on
plot(r2,Mass2/Mass2(end))

grid on;
set(gca,'FontSize',15);
set(gca,'XScale','log')
set(gca,'YScale','log')
