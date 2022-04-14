function [AvRProj,AvR2Proj] = getSersicMoments(nSersic,rHL)
% For derivation of b_n, we use an approximate analytic form from Wikipedia

bnSersic=(-1/3)+(-2194697/30690717750).*nSersic.^(-4)+(131/1148175).* ...
    nSersic.^(-3)+(46/25515).*nSersic.^(-2)+(4/405).*nSersic.^(-1)+2.* ...
    nSersic;

AvRProj = rHL*(bnSersic^(-nSersic))*gamma(3*nSersic)/gamma(2*nSersic);
AvR2Proj = rHL^2*(bnSersic^(-2*nSersic))*gamma(4*nSersic)/gamma(2*nSersic);

end

