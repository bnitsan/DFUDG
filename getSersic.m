function [r,Density,Mass,Sigma,Phi] = getSersic(nSersic,M0Sersic,r0)
% we work in 1e5 Msol/kpc/Gyr units. A small galaxy would have then
% M0Sersic~1e3-1e4, r0~2 (possibly nSersic~0.6)
% For derivation of b_n, we use an approximate analytic form from Wikipedia

G=0.449; % gravitational constant

bnSersic=(-1/3)+(-2194697/30690717750).*nSersic.^(-4)+(131/1148175).* ...
    nSersic.^(-3)+(46/25515).*nSersic.^(-2)+(4/405).*nSersic.^(-1)+2.* ...
    nSersic;

r0Sersic=r0;

rr=0:0.01:22*r0Sersic;
Density = zeros(length(rr),1);
for i=1:length(rr)
    r=rr(i);
    IntRho=@(rp) bnSersic.*exp(1).^((-1).*bnSersic.*(r0Sersic.^(-1).*rp).^( ...
    nSersic.^(-1))).*nSersic.^(-1).*pi.^(-1).*r0Sersic.^(-1).*( ...
    r0Sersic.^(-1).*rp).^((-1)+nSersic.^(-1)).*((-1).*r.^2+rp.^2).^( ...
    -1/2);
    Density(i)=integral(IntRho,r,inf);
end

Mass(1)=0;
%disp(size(rr))
%disp(size(Density))
for i=2:length(rr)
    Mass(i)=trapz(rr(1:i)',4*pi*rr(1:i)'.^2.*Density(1:i));
end

MInf=Mass(end);
Mass = Mass*M0Sersic/MInf;
Density = Density*M0Sersic/MInf;

Phi = zeros(length(rr),1);
intPhi=-G*Mass./(rr.^2);
for i=1:length(rr)-1
    Phi(i)=trapz(rr(i:end),intPhi(i:end));
end
Phi(1) = Phi(2);

intsigma2=G*Density'.*Mass./(rr.^2);
intsigma2(1)=0;

sigma2 = zeros(length(rr),1);
for i=1:length(rr)-1
    sigma2(i)=(1/Density(i))*trapz(rr(i:end),intsigma2(i:end));
end
Sigma = sqrt(sigma2);
Mass = Mass';


%r0Sersic=r0;
r=rr';%*r0Sersic;

end

