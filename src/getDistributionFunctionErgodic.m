function [rr,VFStructsF] = getDistributionFunctionErgodic(r,Density,Mass,Sigma,Phi)
Psi=-Phi;
Vcirc = sqrt(0.449*Mass./r);

dDensitydPsi=diff(Density)./diff(Psi);
%d2DensitydPsi2=diff(dDensitydPsi)./diff(Psi);
dDensitydPsi(1)=1*dDensitydPsi(2);
rInt=r(1:end-1)+diff(r)/2;
PsiInt=Psi(1:end-1)+diff(Psi)/2;

d2DensitydPsi2=diff(dDensitydPsi)./diff(PsiInt);
PsiInt2=PsiInt(1:end-1)+diff(PsiInt)/2;
rInt2=rInt(1:end-1)+diff(rInt)/2;

d3DensitydPsi3=diff(d2DensitydPsi2)./diff(PsiInt2);
PsiInt3=PsiInt2(1:end-1)+diff(PsiInt2)/2;
rInt3=rInt2(1:end-1)+diff(rInt2)/2;

dDensitydPsi   = flip(dDensitydPsi);
PsiInt         = flip(PsiInt);
rInt           = flip(rInt);
d2DensitydPsi2 = flip(d2DensitydPsi2);
PsiInt2        = flip(PsiInt2);
rInt2          = flip(rInt2);
d3DensitydPsi3 = flip(d3DensitydPsi3);
PsiInt3        = flip(PsiInt3);
rInt3          = flip(rInt3);
flipDensity    = flip(Density);
flipSigma      = flip(Sigma);

Eps=PsiInt3;
fEps3= Eps*0;
for jj=2:length(Eps)
    integrand=(2/sqrt(8)/pi^2)*sqrt(-Eps+PsiInt3(jj)).*d3DensitydPsi3;
    fEps3(jj) = trapz(Eps(1:jj),integrand(1:jj));
end

%vesc = zeros(length(rInt3),1);
coarseFactor=10;
VFStructs = cell(round(length(rInt3)/coarseFactor-1),1);
rArray    = zeros(round(length(rInt3)/coarseFactor-1),1);
for ii=1:10:length(rInt3)
    [v,f]=getfvatr(ii,PsiInt3,Eps,fEps3);
    structVF.v = v;
    structVF.f = f;  
    VFStructs{(ii-1)/coarseFactor+1} = structVF;
    rArray((ii-1)/coarseFactor+1)    = rInt3(ii);
end

rr = flip(rArray);
VFStructsF = flip(VFStructs);

end

function [v,f]=getfvatr(ii,PsiInt3,Eps,fEps3)
    PsiR=PsiInt3(ii);
    vrange=sqrt(2*(PsiR-Eps));
    vint = flip(real(vrange(1:ii)));
    fint = flip(fEps3(1:ii));
    vesc(ii) = vint(end);
    
    v=vint;
    f=fint;
end