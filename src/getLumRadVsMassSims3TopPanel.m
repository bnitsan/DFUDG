%%
function getLumRadVsMassSims3TopPanel(MassBins,sims,titleFig,nS,rH,CutRadius,saveFlag,saveName,Vertical)

[DatAvgMass,DatAvgRProj,DatAvgR2Proj,DatNumPerBin]=getDataMoments2(MassBins);


SimsNumPerBin  = zeros(1,length(MassBins)-1);
SimsNum2PerBin = SimsNumPerBin;
SimsAvgRProj   = SimsNumPerBin;
SimsAvgR2Proj  = SimsNumPerBin;
SimsAvgMass    = SimsNumPerBin;

RAvgs  = zeros(length(sims)*3,length(MassBins)-1);
R2Avgs = zeros(length(sims)*3,length(MassBins)-1);
%NAvgs  = zeros(length(sims),length(MassBins)-1);

for ii=1:length(sims)
    [mom{1},mom{2},mom{3}] = getSingleSimMoments(sims{ii},MassBins,CutRadius);
    for jj=1:3    
        SimsNumPerBin = SimsNumPerBin + mom{jj}.num;
        SimsNum2PerBin = SimsNum2PerBin + (mom{jj}.num).^2;
        SimsAvgMass   = (SimsAvgMass.*  (SimsNumPerBin-mom{jj}.num)   + (mom{jj}.m).*mom{jj}.num)./SimsNumPerBin;
        SimsAvgRProj  = (SimsAvgRProj.* (SimsNumPerBin-mom{jj}.num)  + (mom{jj}.r).*mom{jj}.num)./SimsNumPerBin;   
        SimsAvgR2Proj = (SimsAvgR2Proj.*(SimsNumPerBin-mom{jj}.num) + mom{jj}.r2.*mom{jj}.num)./SimsNumPerBin;   

        SimsAvgMass(isnan(SimsAvgMass))=0;
        SimsAvgRProj(isnan(SimsAvgRProj))=0;
        SimsAvgR2Proj(isnan(SimsAvgR2Proj))=0;
    end
    
    RAvgs(3*ii-2,:) = mom{1}.r;
    RAvgs(3*ii-1,:) = mom{2}.r;
    RAvgs(3*ii,:)   = mom{3}.r;
    R2Avgs(3*ii-2,:) = mom{1}.r2;
    R2Avgs(3*ii-1,:) = mom{2}.r2;
    R2Avgs(3*ii,:)   = mom{3}.r2;

end




avR = zeros((length(MassBins)-1),1); upR = avR; downR = avR;
avR2 = zeros((length(MassBins)-1),1); upR2 = avR2; downR2 = avR2;
for ii=1:(length(MassBins)-1)
    [avR(ii),upR(ii),downR(ii)]=getMeanUpDownRad(RAvgs(RAvgs(:,ii)>0,ii));
    [avR2(ii),upR2(ii),downR2(ii)]=getMeanUpDownRad(R2Avgs(R2Avgs(:,ii)>0,ii));
end

% 
%STDr = sqrt(SimsAvgR2Proj-SimsAvgRProj.^2);
STDN = sqrt(SimsNum2PerBin/length(sims)/3-(SimsNumPerBin/length(sims)/3).^2);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get radial initial conditions and mass function initial conditions
% (matching data in case of observed GC mass function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AvRProjIC,AvR2ProjIC] = getSersicMoments(nS,rH);
if sims{1}.MeanInitial == -1
    [AvMassIC,AvNGCIC]=getNGCICObs(MassBins);
elseif sims{1}.MeanInitial == -3
    [AvMassIC,AvNGCIC]=getNGCICAOObs(MassBins);    
elseif sims{1}.MeanInitial == -4
    %[AvMassIC,AvNGCIC]=getNGCICAOObs(MassBins); 
    % fly exception for now in many faint case
    [AvMassIC,AvNGCIC]=getNGCICObs(MassBins);
else
    [AvMassIC,AvNGCIC]=getNGCIC(MassBins,log(sims{1}.MeanInitial),sims{1}.SigmaInitial,sims{1}.N);
end
AvNGCIC(AvNGCIC<0.01)=-1;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% plot figure %%%%%%%%%%%%%
% relevant sim variables: 
% SimsAvgMass,avR,avR-downR,upR-avR
% SimsAvgMass,SimsNumPerBin/length(sims),min(STDN,SimsNumPerBin/length(sims)-0.1),STDN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure('Renderer', 'painters', 'Position', [1200 100 900 400])

figure('Renderer', 'opengl', 'Position', [1200 100 450 450])

scatter(DatAvgMass,DatAvgRProj,120,'filled');
hold on;
errorbar(SimsAvgMass,avR,avR-downR,upR-avR,'.','CapSize',18,'MarkerSize',12,'Color','blue')
plot(0:200,(0:200)*0+AvRProjIC,'linewidth',2)

%axis([0.0 max(max(SimsAvgMass),max(DatAvgMass))*1.1 0 max(max(max(upR),max(DatAvgRProj)),AvRProjIC(1))*1.1])
axis([0.0 27 0.09 5])
xlabel('M bins [$10^5~M_\odot$]','interpreter','latex');
ylabel('$\left<r_\perp\right>$ [kpc]','interpreter','latex')
grid on;
set(gca,'FontSize',15);
set(gca,'YScale','log')

fill([-2 100 100 -2],[1e-4 1e-4 0.5 0.5],'k','FaceAlpha',0.1)
legend('Data','Simulation','Initial condition','Location','NorthEast')
title(titleFig,'interpreter','latex')


if saveFlag
    saveas(gcf,strcat('Plots/',saveName,'TC.eps'),'epsc');
end

end
%%
function [SimMom1,SimMom2,SimMom3] = getSingleSimMoments(sim,MassBins,CutRadius)
mtab = sim.mtab;
XXX = sim.XXX;

m = mtab(end,:);
%binNum=length(MassBins)-1;
%BinDat = ReadDataIntoMassBins(MassBins);

N=sim.N;
i=1:N;
XX = []; YY = []; ZZ = [];
%VX = []; VY = []; VZ = [];
XX(:,i)=XXX(:,i);
YY(:,i)=XXX(:,i+N);
ZZ(:,i)=XXX(:,i+2*N);

%VX(:,i)= XXX(:,i+3*N);
%VY(:,i)= XXX(:,i+4*N);
%VZ(:,i)= XXX(:,i+5*N);

% get simulation result
XM=XX(end,:); YM=YY(end,:); ZM=ZZ(end,:);

rM1=(XM.^2+YM.^2).^0.5;
rM2=(ZM.^2+YM.^2).^0.5;
rM3=(XM.^2+ZM.^2).^0.5;

cond1=rM1<=CutRadius;
cond2=rM2<=CutRadius;
cond3=rM3<=CutRadius;
rM1 = rM1(cond1);
rM2 = rM2(cond2);
rM3 = rM3(cond3);
msub1=m(cond1);
msub2=m(cond2);
msub3=m(cond3);

[mAvg1,rAvg1,r2Avg1,NumPerBinSim1]=getSingleDirectionMoments(msub1,rM1,MassBins);
[mAvg2,rAvg2,r2Avg2,NumPerBinSim2]=getSingleDirectionMoments(msub2,rM2,MassBins);
[mAvg3,rAvg3,r2Avg3,NumPerBinSim3]=getSingleDirectionMoments(msub3,rM3,MassBins);

SimMom1.m = mAvg1;
SimMom1.r = rAvg1;
SimMom1.r2 = r2Avg1;
SimMom1.num = NumPerBinSim1;

SimMom2.m = mAvg2;
SimMom2.r = rAvg2;
SimMom2.r2 = r2Avg2;
SimMom2.num = NumPerBinSim2;

SimMom3.m = mAvg3;
SimMom3.r = rAvg3;
SimMom3.r2 = r2Avg3;
SimMom3.num = NumPerBinSim3;

end

function [mAvg, rAvg1,r2Avg1,NumPerBinSim]=getSingleDirectionMoments(msub1,rM1,MassBins)
MoreMassiveGCs=[];
ind1=find(msub1>MassBins(end));

if ~isempty(ind1)
    MoreMassiveGCs=msub1(ind1);
    disp(strcat('Larger than bin ',num2str(msub1(ind1))))
    msub1(ind1) = [];
    rM1(ind1)   = [];
end

BinSim={};
for ii=1:(length(MassBins)-1)
    ind1=find(msub1>MassBins(end-ii));
    if ~isempty(ind1)
        BinSim{(length(MassBins))-ii}.r1=rM1(ind1);
        BinSim{(length(MassBins))-ii}.mass = msub1(ind1);
        msub1(ind1)  = [];
        rM1(ind1)   = [];
    else
        BinSim{(length(MassBins))-ii}=[1000];
    end
end

mAvg = zeros(1,length(MassBins)-1);
rAvg1 = zeros(1,length(MassBins)-1);
r2Avg1 = zeros(1,length(MassBins)-1);
NumPerBinSim = zeros(1,length(MassBins)-1);

for ii = 1:length(BinSim)
    if isa(BinSim{ii},'double')
        mAvg(ii) = -5;
        rAvg1(ii) = -5;
        r2Avg1(ii) = -5;
        %sizePointSim(ii) = 5;
        NumPerBinSim(ii) = 0;
    else
        mAvg(ii) = mean(BinSim{ii}.mass);
        rAvg1(ii) = mean(BinSim{ii}.r1);
        r2Avg1(ii) = mean((BinSim{ii}.r1).^2);
        %sizePointSim(ii) = 50*log(length(BinSim{ii}.r1)/0.5)^1.5;
        NumPerBinSim(ii) = length(BinSim{ii}.r1);
    end
end
end

function [av,up,down]=getMeanUpDownRad(rAv)

av = mean(rAv);

if length(rAv)>1
    sortedr=sort(rAv);
    la=length(sortedr);
    laD = round(la*0.32/2+0.5);
    laU = round(la*(1-0.32/2)+0.5);

    up = sortedr(laU);
    down = sortedr(laD);
elseif length(rAv)==1
    up = av*2;
    down = av/2;
else
    av=-1;
    up=av;
    down=av;
end
end

function [AvgMass,AvgRProj,AvgR2Proj,NumPerBin]=getDataMoments2(MassBins)

BinDat = ReadDataIntoMassBins(MassBins);

binNum  = length(BinDat);
AvgMass = zeros(binNum,1);
%ErrMass = zeros(binNum,1);
AvgRProj = zeros(binNum,1);
AvgR2Proj = zeros(binNum,1);
sizePoint = zeros(binNum,1);
NumPerBin = zeros(binNum,1);
for ii = 1:binNum
    if isempty(BinDat{ii})
        disp('yo bin empty')
        AvgMass(ii) = -1;
        %ErrMass(ii) = std(BinDat{ii}(:,1));
        AvgRProj(ii) = -1;
        AvgR2Proj(ii) = -1;
        sizePoint(ii) = -1;
        NumPerBin(ii) = 0;
    else
        AvgMass(ii) = mean(BinDat{ii}(:,2));
        %ErrMass(ii) = std(BinDat{ii}(:,1));
        AvgRProj(ii) = mean(BinDat{ii}(:,1));
        AvgR2Proj(ii) = mean((BinDat{ii}(:,1)).^2);
        sizePoint(ii) = 50*log(length(BinDat{ii})/0.5)^1.5;
        NumPerBin(ii) = length(BinDat{ii}(:,1));
    end
    
end

end

function [mAvg,rAvg1,rAvg2,rAvg3,r2Avg1,r2Avg2,r2Avg3,NumPerBinSim]=getSimMoments3(sim,MassBins)
mtab = sim.mtab;
XXX = sim.XXX;

m = mtab(end,:);
%binNum=length(MassBins)-1;
%BinDat = ReadDataIntoMassBins(MassBins);

N=sim.N;
i=1:N;
XX = []; YY = []; ZZ = [];
%VX = []; VY = []; VZ = [];
XX(:,i)=XXX(:,i);
YY(:,i)=XXX(:,i+N);
ZZ(:,i)=XXX(:,i+2*N);

%VX(:,i)= XXX(:,i+3*N);
%VY(:,i)= XXX(:,i+4*N);
%VZ(:,i)= XXX(:,i+5*N);

% get simulation result
XM=XX(end,:); YM=YY(end,:); ZM=ZZ(end,:);

rM1=(XM.^2+YM.^2).^0.5;
rM2=(ZM.^2+YM.^2).^0.5;
rM3=(XM.^2+ZM.^2).^0.5;

MoreMassiveGCs=[];
msub=m;
ind=find(msub>MassBins(end));

if ~isempty(ind)
    MoreMassiveGCs=msub(ind);
    disp(strcat('Larger than bin ',num2str(msub(ind))))
    msub(ind) = [];
    rM1(ind)   = [];
    rM2(ind)   = [];
    rM3(ind)   = [];    
end

BinSim={};
%effLenBin = length(MassBins)-1;
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
        %effLenBin = (length(MassBins))-2;
    end
end

mAvg = zeros(1,length(MassBins)-1);
rAvg1 = zeros(1,length(MassBins)-1);
rAvg2 = zeros(1,length(MassBins)-1);
rAvg3 = zeros(1,length(MassBins)-1);
r2Avg1 = zeros(1,length(MassBins)-1);
r2Avg2 = zeros(1,length(MassBins)-1);
r2Avg3 = zeros(1,length(MassBins)-1);

for ii = 1:length(BinSim)
    if isa(BinSim{ii},'double')
        mAvg(ii) = -5;
        rAvg(ii) = -5;
        r2Avg(ii) = -5;
        sizePointSim(ii) = 5;
        NumPerBinSim(ii) = 0;
    else
        mAvg(ii) = mean(BinSim{ii}.mass);
        rAvg1(ii) = mean(BinSim{ii}.r1);
        rAvg2(ii) = mean(BinSim{ii}.r2);
        rAvg3(ii) = mean(BinSim{ii}.r3);
        r2Avg1(ii) = mean((BinSim{ii}.r1).^2);
        r2Avg2(ii) = mean((BinSim{ii}.r2).^2);
        r2Avg3(ii) = mean((BinSim{ii}.r3).^2);
        %r2Avg(ii) = mean((BinSim{ii}.r1).^2+(BinSim{ii}.r2).^2+(BinSim{ii}.r3).^2)/3;
        %sizePointSim(ii) = 50*log(length(BinSim{ii}.r1)/0.5)^1.5;
        NumPerBinSim(ii) = length(BinSim{ii}.r1);
    end
end

end

function [mAvg,rAvg1,rAvg2,rAvg3,NumPerBinSim]=getSimMoments2(sim,MassBins)
mtab = sim.mtab;
XXX = sim.XXX;

m = mtab(end,:);
%binNum=length(MassBins)-1;
%BinDat = ReadDataIntoMassBins(MassBins);

N=sim.N;
i=1:N;
XX = []; YY = []; ZZ = [];
%VX = []; VY = []; VZ = [];
XX(:,i)=XXX(:,i);
YY(:,i)=XXX(:,i+N);
ZZ(:,i)=XXX(:,i+2*N);

%VX(:,i)= XXX(:,i+3*N);
%VY(:,i)= XXX(:,i+4*N);
%VZ(:,i)= XXX(:,i+5*N);

% get simulation result
XM=XX(end,:); YM=YY(end,:); ZM=ZZ(end,:);

rM1=(XM.^2+YM.^2).^0.5;
rM2=(ZM.^2+YM.^2).^0.5;
rM3=(XM.^2+ZM.^2).^0.5;

MoreMassiveGCs=[];
msub=m;
ind=find(msub>MassBins(end));

if ~isempty(ind)
    MoreMassiveGCs=msub(ind);
    disp(strcat('Larger than bin ',num2str(msub(ind))))
    msub(ind) = [];
    rM1(ind)   = [];
    rM2(ind)   = [];
    rM3(ind)   = [];    
end

BinSim={};
%effLenBin = length(MassBins)-1;
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
        %effLenBin = (length(MassBins))-2;
    end
end

mAvg = zeros(1,length(MassBins)-1);
rAvg1 = zeros(1,length(MassBins)-1);
rAvg2 = zeros(1,length(MassBins)-1);
rAvg3 = zeros(1,length(MassBins)-1);

for ii = 1:length(BinSim)
    if isa(BinSim{ii},'double')
        mAvg(ii) = -5;
        rAvg(ii) = -5;
        r2Avg(ii) = -5;
        sizePointSim(ii) = 5;
        NumPerBinSim(ii) = 0;
    else
        mAvg(ii) = mean(BinSim{ii}.mass);
        rAvg1(ii) = mean(BinSim{ii}.r1);
        rAvg2(ii) = mean(BinSim{ii}.r2);
        rAvg3(ii) = mean(BinSim{ii}.r3);
        %r2Avg(ii) = mean((BinSim{ii}.r1).^2+(BinSim{ii}.r2).^2+(BinSim{ii}.r3).^2)/3;
        %sizePointSim(ii) = 50*log(length(BinSim{ii}.r1)/0.5)^1.5;
        NumPerBinSim(ii) = length(BinSim{ii}.r1);
    end
end

end