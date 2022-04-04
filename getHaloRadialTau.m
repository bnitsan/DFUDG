function [rS0,TausS0,rB,TausB,rC,TausC] = getHaloRadialTau()

SersicMass = 10*1.2e3; SersicRadius = 1.9; nSersicStars=0.61;
%[rS,DensityS,MassS,SigmaS,PhiS] = getSersic(nSersicStars,SersicMass,SersicRadius);
%VcircS = sqrt(0.449*MassS./rS);
[rS0,DensityS0,MassS0,SigmaS0,PhiS0] = getSersic(nSersicStars,0.1*SersicMass,SersicRadius);
VcircS0 = sqrt(0.449*MassS0./rS0);


rsNFW=6; cNFW = 6;
%rsNFW=5; cNFW = 9;
[rC,DensityC,MassC,SigmaC,PhiC] = getNFW(rsNFW,cNFW);
VcircC = sqrt(0.449*MassC./rC);

rsBur = 2; cBur = 15;
[rB,DensityB,MassB,SigmaB,PhiB] = getBurkert(rsBur,cBur);
VcircB = sqrt(0.449*MassB./rB);


%XvS  = VcircS./SigmaS/sqrt(2);
XvS0  = VcircS0./SigmaS0/sqrt(2);
XvC  = VcircC./SigmaC/sqrt(2);
XvB  = VcircB./SigmaB/sqrt(2);

m=5;
TausB        = 1.*0.789.*VcircB.^3./(m.*DensityB.*...
    (erf(XvB)-2*XvB.*exp(-XvB.^2)/sqrt(pi)).*...
    log(1+(4.45*rB.*VcircB.^2./m).^2));
TausS0        = 1.*0.789.*VcircS0.^3./(m.*DensityS0.*...
    (erf(XvS0)-2*XvS0.*exp(-XvS0.^2)/sqrt(pi)).*...
    log(1+(4.45*rS0.*VcircS0.^2./m).^2));
TausC        = 1.*0.789.*VcircC.^3./(m.*DensityC.*...
    (erf(XvC)-2*XvC.*exp(-XvC.^2)/sqrt(pi)).*...
    log(1+(2.227*0.5.*SigmaC.^2./m).^2));

end

