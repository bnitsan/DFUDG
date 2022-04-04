function [r,projCDF] = getSersicProjCDF(nSersic,r0)
% we work in 1e5 Msol/kpc/Gyr units. A small galaxy would have then
% M0Sersic~1e3-1e4, r0~2 (possibly nSersic~0.6)
% For derivation of b_n, we use an approximate analytic form from Wikipedia

G=0.449; % gravitational constant

bnSersic=(-1/3)+(-2194697/30690717750).*nSersic.^(-4)+(131/1148175).* ...
    nSersic.^(-3)+(46/25515).*nSersic.^(-2)+(4/405).*nSersic.^(-1)+2.* ...
    nSersic; 

r0Sersic=r0;

rr=0:0.005:22*r0Sersic;

SigmaDensity = exp(-bnSersic*(rr/r0Sersic).^(1/nSersic));
projCDFArb = zeros(length(rr),1);
for ii = 2:length(rr)
    projCDFArb(ii) = projCDFArb(ii-1)+ trapz(rr((ii-1):ii),rr((ii-1):ii).*SigmaDensity((ii-1):ii));
end

projCDF = projCDFArb./(projCDFArb(end));

r=rr';

end