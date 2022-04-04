function [BinDat] = ReadDataIntoMassBins(MassBins)

MassUnits=1e5;
MassToLight = 1.6;
Distance = 26.5e3; % kpc
dat = ReadData();

%tab = [xid id RA Dec f606tot f475tot vel vel_err specconf f606err f475err];
% distance
RACenter  = 226.3342;
DecCenter = 1.8128;
rperp     = Distance*(pi/180)*sqrt((dat(:,3)-RACenter).^2.*(cos((pi/180/2)*(dat(:,4)+DecCenter))).^2+(dat(:,4)-DecCenter).^2);%mass =

% luminosity/mass
MX606     = dat(:,5) - 5*log(Distance*1e3/10)/log(10);
L606      = 10.^((MX606-4.8)/(-2.5));
MGC       = MassToLight*L606/MassUnits;

% bins
ind=find(MGC>MassBins(end));
if ~isempty(ind)
    [rperp(ind) MGC(ind)];
    rperp(ind) = [];
    MGC(ind)   = [];
end

BinDat={};
for ii=1:(length(MassBins)-1)
    %disp(strcat(num2str(MassBins(end-ii+1)),'and',num2str(MassBins(end-ii))))
    ind=find(MGC>MassBins(end-ii));
    if ~isempty(ind)
        BinDat{(length(MassBins))-ii}=[rperp(ind) MGC(ind)];
        rperp(ind) = [];
        MGC(ind)   = [];
    end
end


end

