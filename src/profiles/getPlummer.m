function [r,Density,Mass,Sigma,Phi] = getPlummer(MPlummer,rPlummer)
% we work in 1e5 Msol/kpc/Gyr units. A small galaxy would have then
% MPlummer~1e3-1e4, rPlummer~2

G=0.449; % gravitational constant

r=0:0.005:22*2;
a=rPlummer;
Mtot=MPlummer;
Density=3*(Mtot)*(1/4/pi/a^3)*(1+r.^2/a^2).^(-5/2);
Mass=Mtot*r.^3./(r.^2+a^2).^(3/2);
Phi =-G*Mtot./sqrt(r.^2+a^2);
Sigma = sqrt(G*Mtot*(1/6)./sqrt(r.^2+a^2));

end

