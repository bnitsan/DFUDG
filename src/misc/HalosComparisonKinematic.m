clear
%%
SersicMass = 1.2e3; SersicRadius = 1.9; nSersicStars=0.61;
%[rS,DensityS,MassS,SigmaS,PhiS] = getSersic(nSersicStars,SersicMass,SersicRadius);
%VcircS = sqrt(0.449*MassS./rS);
[rS0,DensityS0,MassS0,SigmaS0,PhiS0] = getSersic(nSersicStars,SersicMass,SersicRadius);
VcircS0 = sqrt(0.449*MassS0./rS0);


rsNFW=6; cNFW = 6;%cNFW = 6;
%rsNFW=5; cNFW = 9;
[rC,DensityC,MassC,SigmaC,PhiC] = getNFW(rsNFW,cNFW);
VcircC = sqrt(0.449*MassC./rC);

rsBur = 2; cBur = 15;
[rB,DensityB,MassB,SigmaB,PhiB] = getBurkert(rsBur,cBur);
VcircB = sqrt(0.449*MassB./rB);

%%
%XvS  = VcircS./SigmaS/sqrt(2);
XvS0  = VcircS0./SigmaS0/sqrt(2);
XvC  = VcircC./SigmaC/sqrt(2);
XvB  = VcircB./SigmaB/sqrt(2);

m=5;
TausB        = 1.*0.789.*VcircB.^3./(m.*DensityB.*...
    (erf(XvB)-2*XvB.*exp(-XvB.^2)/sqrt(pi)).*...
    log(1+(4.45*rB.*VcircB.^2./m).^2));

%SmoothHeaviside = 1./(1+exp(-(rB-2)/0.1));
%TausB        = 1.*0.789.*VcircB.^3./(m.*DensityB.*...
%    (erf(XvB)-2*XvB.*exp(-XvB.^2)/sqrt(pi)).*...
%    (log(1+(4.45*rB.*VcircB.^2./m).^2)+(log(1+(2.227*0.5.*SigmaB.^2./m).^2)-log(1+(4.45*rB.*VcircB.^2./m).^2)).*SmoothHeaviside));

TausS0        = 1.*0.789.*VcircS0.^3./(m.*DensityS0.*...
    (erf(XvS0)-2*XvS0.*exp(-XvS0.^2)/sqrt(pi)).*...
    log(1+(4.45*rS0.*VcircS0.^2./m).^2));
TausC        = 1.*0.789.*VcircC.^3./(m.*DensityC.*...
    (erf(XvC)-2*XvC.*exp(-XvC.^2)/sqrt(pi)).*...
    log(1+(2.227*0.5.*SigmaC.^2./m).^2));

%%
figure
plot(rS0, DensityS0,'--','linewidth',2)
hold on;
%plot(rS, DensityS,'linewidth',2)
plot(rB, DensityB,'-.','linewidth',2)

plot(rC, DensityC,'linewidth',2)

set(gca,'XScale','log')
set(gca,'YScale','log')
grid on;
set(gca,'FontSize',15);
xlabel('r [kpc]','interpreter','latex');
ylabel('Density [$10^5~M_\odot$/kpc$^3$]','interpreter','latex');
axis([0.1 5 1e0 2e3])
legend('Stars','Burkert','NFW')
%saveas(gcf,'Plots/UDGDensity.eps')
saveas(gcf,strcat('Plots/UDGDensity.eps'),'epsc');

%%
figure
plot(rS0, MassS0,'--','linewidth',2)
hold on;
plot(rB, MassB,'-.','linewidth',2)
plot(rC, MassC,'linewidth',2)

dat=ReadDataIntoMassBins([0.1 30]);
datSort = sortrows(dat{1},1);
for ii=2:length(datSort)
    datSort(ii,2) = datSort(ii,2)+datSort(ii-1,2);
end
stairs([0.0001 (5/4)*datSort(:,1)'],[0.0001 datSort(:,2)'],'linewidth',2)

set(gca,'XScale','log')
set(gca,'YScale','log')
grid on;
set(gca,'FontSize',15);
xlabel('r [kpc]','interpreter','latex');
ylabel('M(r) [$10^5~M_\odot$]','interpreter','latex');
axis([0.1 5 1e0 1e4])
legend('Stars','Burkert','NFW','GC mass','Location','NorthWest')
%saveas(gcf,'Plots/UDGmassEnc.pdf')
saveas(gcf,strcat('Plots/UDGmassEnc.eps'),'epsc');

%%
figure
plot(rS0, 0.979*VcircS0,'--','linewidth',2)
hold on;
plot(rB, 0.979*VcircB,'-.','linewidth',2)

plot(rC, 0.979*VcircC,'linewidth',2)
axis([0 4 0 40])
grid on;
set(gca,'FontSize',15);
xlabel('r [kpc]','interpreter','latex');
ylabel('$V_{\rm circ}$ [km/sec]','interpreter','latex');
legend('Stars','Burkert','NFW','Location','NorthWest')
%saveas(gcf,'Plots/UDGVcirc.pdf')
saveas(gcf,strcat('Plots/UDGVcirc.eps'),'epsc');

%%
figure
plot(rS0, 0.979*SigmaS0,'--','linewidth',2)
hold on;
plot(rB, 0.979*SigmaB,'-.','linewidth',2)

plot(rC, 0.979*SigmaC,'linewidth',2)
axis([0 4 0 27])
grid on;
set(gca,'FontSize',15);
xlabel('r [kpc]','interpreter','latex');
ylabel('$\sigma$ [km/sec]','interpreter','latex');
legend('Stars','Burkert','NFW','Location','NorthEast')
%saveas(gcf,'Plots/UDGSigma.pdf')
saveas(gcf,strcat('Plots/UDGSigma.eps'),'epsc');

%%
rStall=0.2;
TauFudgeFactorOverall=1;
TauFudgeFac=((rC<rStall)*10+(rC>=rStall)*1)*TauFudgeFactorOverall;
TauFudgeFacAdd=((rC<rStall)*10+(rC>=rStall)*0)*TauFudgeFactorOverall;

figure
plot(rS0, TausS0,'--','linewidth',2)
hold on;
plot(rB, TausB,'-.','linewidth',2)

%plot(rC, TausC,'linewidth',2)
plot(rC, TausC+TauFudgeFacAdd,'linewidth',2)
%plot(rC, 10*(rC/1.15).^2,'linewidth',2)
set(gca,'YScale','log')
set(gca,'XScale','log')

axis([0.2 5 0.3 67])
grid on;
set(gca,'FontSize',15);
xlabel('r [kpc]','interpreter','latex');
ylabel('$\tau_{\rm DF}(r,V_{\rm circ}(r))$ [Gyr]','interpreter','latex');
legend('Stars','Burkert','NFW','Location','SouthEast')
%saveas(gcf,'Plots/UDGTauDF.eps')
saveas(gcf,strcat('Plots/UDGTauDF.eps'),'epsc');

%%
figure
plot(rS0, PhiS0,'--','linewidth',2)
hold on;
plot(rB, PhiB,'-.','linewidth',2)

plot(rC, PhiC,'linewidth',2)
axis([0 4 -4e3 0])
grid on;
set(gca,'FontSize',15);
xlabel('$r$ [kpc]','interpreter','latex');
ylabel('$\Phi$','interpreter','latex');
legend('Stars','Burkert','NFW','Location','SouthEast')
%saveas(gcf,'Plots/UDGPhiExt.pdf')
saveas(gcf,strcat('Plots/UDGPhiExt.eps'),'epsc');

%%
nSersicGC = 0.61;
SersicRadiusGC = 1.9;
SersicMassGC = 130;
[rGC,DensityGC,MassGC,SigmaGC,PhiGC] = getSersic(nSersicGC,SersicMassGC,SersicRadiusGC);
bnSersicGCs=(-1/3)+(-2194697/30690717750).*nSersicGC.^(-4)+(131/1148175).* ...
    nSersicGC.^(-3)+(46/25515).*nSersicGC.^(-2)+(4/405).*nSersicGC.^(-1)+2.* ...
    nSersicGC;
NormI = SersicMassGC/(2*pi*nSersicGC*SersicRadiusGC^2*(bnSersicGCs^(-2*nSersicGC))*gamma(2*nSersicGC));
Ir = NormI*exp(-bnSersicGCs*(rGC/SersicRadiusGC).^(1/nSersicGC));

%%
nuSigmaB2 = DensityGC.*interp1(rB,SigmaB,rGC).^2;
nuSigmaS02 = DensityGC.*interp1(rS0,SigmaS0,rGC).^2;
nuSigmaC2 = DensityGC.*interp1(rC,SigmaC,rGC).^2;

dnusigma2drpB=[(diff(nuSigmaB2)./diff(rGC))' 0]';
dnusigma2drpS0=[(diff(nuSigmaS02)./diff(rGC))' 0]';
dnusigma2drpC=[(diff(nuSigmaC2)./diff(rGC))' 0]';

rSigmaLOS = rGC;
sLosIntB = zeros(length(rSigmaLOS)-1,1);
sLosIntS0 = zeros(length(rSigmaLOS)-1,1);
sLosIntC = zeros(length(rSigmaLOS)-1,1);
for ii = 1:length(rSigmaLOS)-1
    arrB  =-dnusigma2drpB.*sqrt(rGC.^2-rSigmaLOS(ii)^2);
    arrS0 =-dnusigma2drpS0.*sqrt(rGC.^2-rSigmaLOS(ii)^2);
    arrC  =-dnusigma2drpC.*sqrt(rGC.^2-rSigmaLOS(ii)^2);
    sLosIntB(ii)=trapz(rGC(ii:end),arrB(ii:end));
    sLosIntS0(ii)=trapz(rGC(ii:end),arrS0(ii:end));
    sLosIntC(ii)=trapz(rGC(ii:end),arrC(ii:end));
end
sLosIntB = [sLosIntB' 0]';
sLosIntS0 = [sLosIntS0' 0]';
sLosIntC = [sLosIntC' 0]';
%%
SigmaBLOS=sqrt((2./Ir).*sLosIntB);
SigmaS0LOS=sqrt((2./Ir).*sLosIntS0);
SigmaCLOS=sqrt((2./Ir).*sLosIntC);

figure
plot(rGC, 0.979*SigmaS0LOS,'--','linewidth',2)
hold on
plot(rGC, 0.979*SigmaBLOS,'-.','linewidth',2)
plot(rGC, 0.979*SigmaCLOS,'linewidth',2)

grid on;
set(gca,'FontSize',15);
xlabel('$r_\perp$ [kpc]','interpreter','latex');
ylabel('$\sigma_{\rm LOS}$ [km/sec]','interpreter','latex');
legend('Stars','Stars*10 (mass)','NFW','Location','NorthEast')
axis([0.0 4.5 0 20])

legend('Stars','Burkert','NFW','Location','SouthWest')
%saveas(gcf,'Plots/UDGSigmaLOS.pdf')
saveas(gcf,strcat('Plots/UDGSigmaLOS.eps'),'epsc');
