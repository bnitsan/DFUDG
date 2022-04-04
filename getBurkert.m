function [r,Density,Mass,Sigma,Phi] = getBurkert(rs,c)
% a function that gets a 1) radius - rs, 2) concentration - c
% and returns [radius, density, enclosed mass, radial velocity dispersion, gravitational potential]
% corresponding a Burkert profile.
%
% we work in 1e5 Msol/kpc/Gyr units. A small galaxy would have then
% M0Sersic~1e3-1e4, r0~2 (possibly nSersic~0.6)

G=0.449; % gravitational constant

KmSToKpcGyr = 1/0.979;

rhoc   = 3*(70*KmSToKpcGyr/1000)^2/(8*pi*G);
deltac = (200/3)*c^3/(log(1+c)-c/(1+c));
rhos   = deltac*rhoc;

rr=0:0.005:22*rs;

Density = rhos./(((rr+rs)/rs).*(1+(rr/rs).^2));

Mass = zeros(1,length(rr));
for ii=2:length(rr)
    Mass(ii)=Mass(ii-1)+trapz(rr((ii-1):ii),4*pi*rr((ii-1):ii).^2.*Density((ii-1):ii));
end

rr2=rr(2:end);
Phi(2:length(rr)) = G*rhos*rs^2*(-pi./rr2).*(-rr2.*log(1+1./rr2.^2)+...
    log(1+rr2.^2)+2.*log(rr2+1)-2*atan(rr2)+...
    2*rr2.*acot(rr2)+4*rr2.*acoth(1+2.*rr2));

Phi(1) = Phi(2); % fix singularity

intsigma2=G*Density.*Mass./(rr.^2);
intsigma2(1)=0;

sigma2 = zeros(length(rr),1);
for ii=(length(rr)-1):(-1):1
    sigma2(ii) = sigma2(ii+1)*Density(ii+1)/Density(ii) + (1/Density(ii))*trapz(rr(ii:(ii+1)),intsigma2(ii:(ii+1)));
end
Sigma = sqrt(sigma2);
Sigma = Sigma';
r=rr;

end

