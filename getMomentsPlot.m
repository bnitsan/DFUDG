function h=getMomentsPlot(sim,MassBins)

mtab = sim.mtab;
XXX = sim.XXX;

m = mtab(end,:);

N=sim.N;
i=1:N;
XX = []; YY = []; ZZ = [];

XX(:,i)=XXX(:,i);
YY(:,i)=XXX(:,i+N);
ZZ(:,i)=XXX(:,i+2*N);

% get simulation result
XM=XX(end,:); YM=YY(end,:); ZM=ZZ(end,:);

rM1=(XM.^2+YM.^2).^0.5;
rM2=(ZM.^2+YM.^2).^0.5;
rM3=(XM.^2+ZM.^2).^0.5;

msub=m;
ind=find(msub>MassBins(end));

if ~isempty(ind)
    disp(strcat('Larger than bin ',num2str(msub(ind))))
    msub(ind) = [];
    rM1(ind)   = [];
    rM2(ind)   = [];
    rM3(ind)   = [];    
end

BinSim={};

for ii=1:(length(MassBins)-1)

    ind=find(msub>MassBins(end-ii));
    if ~isempty(ind)
        BinSim{(length(MassBins))-ii}.r1=rM1(ind);
        BinSim{(length(MassBins))-ii}.r2=rM2(ind);
        BinSim{(length(MassBins))-ii}.r3=rM3(ind);
        BinSim{(length(MassBins))-ii}.mass = msub(ind);
        msub(ind)  = [];
        rM1(ind)   = [];
        rM2(ind)   = [];
        rM3(ind)   = [];    
    else
        BinSim{(length(MassBins))-ii}=[1000];
    end
end

for ii = 1:length(BinSim)
    if isa(BinSim{ii},'double')
        mAvg(ii) = -5;
        rAvg(ii) = -5;
        sizePointSim(ii) = 5;
        NumPerBinSim(ii) = 0;
    else
        mAvg(ii) = mean(BinSim{ii}.mass);
        rAvg(ii) = mean(BinSim{ii}.r1+BinSim{ii}.r2+BinSim{ii}.r3)/3;
        sizePointSim(ii) = 50*log(length(BinSim{ii}.r1)/0.5)^1.5;
        NumPerBinSim(ii) = length(BinSim{ii}.r1);
    end
end

h=plotGCMassvsrprojMoment(MassBins);
scatter(rAvg,mAvg,sizePointSim,'filled')
legend('Data','Data-driven model','Simulation')
text(rAvg',mAvg'+1,num2str(NumPerBinSim'),'FontSize',14)

axis([0.1 4 0 MassBins(end)])
end

