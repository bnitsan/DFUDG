function [AvMass,AvNGC]=getNGCIC(MassBins,mu,sigma,Nin)
% first do masses expectation

%mu=log(3);
%sigma=0.3;
fM = @(M) exp(-0.5*(log(M)-mu).^2/sigma^2)./M;
fMM = @(M) M.*exp(-0.5*(log(M)-mu).^2/sigma^2)./M;

AvMass = zeros(length(MassBins)-1,1);
for ii = 1:(length(MassBins)-1)
    AvMass(ii)=integral(fMM,MassBins(ii),MassBins(ii+1))/integral(fM,MassBins(ii),MassBins(ii+1));
end
% then do N_GC

AvNGC = zeros(length(MassBins)-1,1);
for ii = 1:(length(MassBins)-1)
    AvNGC(ii)=integral(fM,MassBins(ii),MassBins(ii+1))/integral(fM,0,inf);
end

AvNGC=AvNGC*Nin;
end