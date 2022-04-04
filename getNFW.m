function [r,Density,Mass,Sigma,Phi] = getNFW(rs,c)
% a function that gets a 1) radius - rs, 2) concentration - c
% and returns [radius, density, enclosed mass, radial velocity dispersion, gravitational potential]
% corresponding a NFW profile.
%
% we work in 1e5 Msol/kpc/Gyr units.

G=0.449; % gravitational constant

KmSToKpcGyr = 1/0.979;

rhoc   = 3*(70*KmSToKpcGyr/1000)^2/(8*pi*G);
deltac = (200/3)*c^3/(log(1+c)-c/(1+c));
rhos   = deltac*rhoc;

rr=0:0.005:22*rs;

Density = rhos./((rr/rs).*(1+rr/rs).^2);
Density(1) = Density(2);

Mass = zeros(1,length(rr));
for ii=2:length(rr)
    Mass(ii)=Mass(ii-1)+trapz(rr((ii-1):ii),4*pi*rr((ii-1):ii).^2.*Density((ii-1):ii));
end

Phih   = -4*pi*G*rhos*rs^2;
Phi(1) = Phih;
Phi(2:length(rr)) = Phih.*(rs./rr(2:end)).*log(1+rr(2:end)/rs);

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

