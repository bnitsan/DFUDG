function [AvMass,NGCarr]=getNGCICAOObs(MassBins)
% first do masses expectation

datM = getGCMassSpecAlmostObserved();

%dat=ReadDataIntoMassBins(MassBins);

NGCarr = zeros(length(MassBins)-1,1)-1;
AvMass = zeros(length(MassBins)-1,1);

for ii = 2:length(MassBins)
    datMtemp = datM(datM<(MassBins(ii)));
    datMF    = datMtemp(datMtemp>(MassBins(ii-1)));
    if ~isempty(datMF)
        NGCarr(ii-1) = length(datMF);
        AvMass(ii-1) = mean(datMF);
    end
end
%disp(num2str(datMF))
%for ii = 1:length(dat)
%    NGCarr(ii) = size(dat{ii},1);
%    AvMass(ii) = mean(dat{ii}(:,2));
%end
end