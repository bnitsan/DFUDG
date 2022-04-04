function [AvMass,NGCarr]=getNGCICObs(MassBins)
% first do masses expectation

dat=ReadDataIntoMassBins(MassBins);
NGCarr = zeros(length(MassBins)-1,1)-1;
AvMass = zeros(length(MassBins)-1,1);
for ii = 1:length(dat)
    if isempty(dat{ii})
        NGCarr(ii) = -1;
        AvMass(ii) = -1;
    else
        NGCarr(ii) = size(dat{ii},1);
        AvMass(ii) = mean(dat{ii}(:,2));
    end
end
end