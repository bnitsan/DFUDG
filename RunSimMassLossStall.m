function SimStruct = RunSimMassLossStall(HaloParameters,SersicRadiusGCs,nSersicGCs,Nin,SimTime,GnFrac,TauFudgeFactorOverall,IsoVsCircFlag,SofteningRadKpc,massScaleIn,MassNormMean,MassNormVar,MergeRadius,TidalRadius,toleranceFactor,dtFactor,MGCMhaloRestrictFraction,TotalLossFraction)
% a single simulation to use in an ensemble of simulations
global G Gn m MassLossPerStep LogOption massscale TauFudgeFac critRadius MergeList sigmaS rS densityS massS phiS epsilon maxrc rTidal rmerged realMergeTidal

rSersicGCs = SersicRadiusGCs;
MSersicGCs = 1e3; % this is irrelevant, really

% define halo properties from input config
if HaloParameters.type == 1
    [rS,densityS,massS,sigmaS,phiS]=getSersic(HaloParameters.nSersic,HaloParameters.MSersic,HaloParameters.rSersic);
    LogOption = 1;
elseif HaloParameters.type == 2
    [rS,densityS,massS,sigmaS,phiS]=getNFW(HaloParameters.NFWRs,HaloParameters.NFWc);
    LogOption = 2;
elseif HaloParameters.type == 3
    [rS,densityS,massS,sigmaS,phiS]=getBurkert(HaloParameters.BurRs,HaloParameters.Burc);
    LogOption = 1;
end
[rSGCs,densitySGCs,massSGCs,~,~]=getSersic(nSersicGCs,MSersicGCs,rSersicGCs);

G=0.449; % gravitational constant
N=Nin; % number of objects
Gn = G*GnFrac; % the gravitational constant between point masses; allows us to turn off nbody
IsotropicVsCircular=IsoVsCircFlag; % assume 0 = isotropic, 1 = circular
epsilon = SofteningRadKpc; % softening radius for Plummer model of GCs

% initialize a stalling factor of dynamical friction
TauInStall = 50;
rStall=HaloParameters.rStall;
TauFudgeFac=((rS<rStall)*TauInStall+(rS>=rStall)*0)*TauFudgeFactorOverall; % this factor gets added to tau, no multiplied

% initalize parameters of mass loss
TimeStepMassLoss = 0.1;
MassLossFractionPerStep = TimeStepMassLoss*TotalLossFraction;
IntegratedLoss = TotalLossFraction*SimTime;

% initialize GC mass function and mass loss
if length(massScaleIn)>1 % if massScaleIn represents already an input GC mass list 
    massscale=5; % arbitrary
    m = massScaleIn;
    m = m/(1-IntegratedLoss);
    MassLossPerStep = MassLossFractionPerStep*m;
else                     % if a Gaussian distribution is assumed
    massscale=massScaleIn;
    m = SampleMassFunction(MassNormMean,MassNormVar,N);
    m = m/(1-IntegratedLoss);
    MassLossPerStep = MassLossFractionPerStep*m;
end

% subtract extra GC mass from observed star mass in the case of Stars-only halo
if HaloParameters.type == 1 
    MLRatio  = (massS(end)-sum(IntegratedLoss*m))/massS(end);
    massS    = MLRatio*massS;
    densityS = MLRatio*densityS;
    sigmaS   = MLRatio*sigmaS;
    phiS     = MLRatio*phiS;    
end

critRadius=MergeRadius; % radius below which GCs may merge
rTidal=TidalRadius;     % radius above which GCs assumed to be tidally disrupted

Rchar = 1*rSersicGCs; % charcteristic radius, used in timescale estimate of integration
rmerged = TidalRadius-0.1;   % where to put GCs after merger





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE RANDOM PHASE-SPACE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% POSITIONS
% azimuthal angles
phi      = rand(N,1)*2*pi;
costheta = rand(N,1)*2-1;

% radial sampling
SersicCDF = massSGCs/massSGCs(end);
iCDFEnd=find(SersicCDF>=0.999999999,1,'first');
rs=interp1(SersicCDF(1:iCDFEnd),rSGCs(1:iCDFEnd),rand(N,1));

% Cartesian coordinates of the objects
X = rs.*cos(phi).*sqrt(1-costheta.^2);         
Y = rs.*sin(phi).*sqrt(1-costheta.^2);        
Z = rs.*costheta;               

rc0=sqrt(X.^2+Y.^2+Z.^2);

r=[X Y Z];                  % position vectors of all masses

% VELOCITIES
if IsotropicVsCircular==0   % 0 - isotropic distribution
    GCDensityInterp1 = interp1(rSGCs,densitySGCs,rS);
    GCDensityInterp1(isnan(GCDensityInterp1))=0;
    [rr,VFStructsF] = getDistributionFunctionErgodic(rS,GCDensityInterp1,massS,sigmaS,phiS);
    
    vRand = zeros(length(rc0),1);
    for kk = 1:length(rc0)
        ix=find(rr>rc0(kk),1,'first');
        totCDF = trapz(VFStructsF{ix}.v,4*pi*(VFStructsF{ix}.v).^2.*VFStructsF{ix}.f);
        xx=VFStructsF{ix}.v;
        yy=VFStructsF{ix}.f;
        CDFLocal = zeros(length(xx),1);

        for jj = 2:length(xx)
            CDFLocal(jj) = CDFLocal(jj-1) + trapz(xx((jj-1):jj),4*pi*xx((jj-1):jj).^2.*yy((jj-1):jj))/totCDF;
        end
        
        randUnif=rand(1);
        iRand=find(CDFLocal>randUnif,1,'first');
        vRand(kk) = xx(iRand);
    end

    phiV      = rand(N,1)*2*pi;
    costhetaV = rand(N,1)*2-1;

    vx = vRand.*cos(phiV).*sqrt(1-costhetaV.^2);         % Cartesian coordinates of the objects
    vy = vRand.*sin(phiV).*sqrt(1-costhetaV.^2);        
    vz = vRand.*costhetaV;               

else % circular velocity distribution
    vcircr   = sqrt(G*interp1(rS,massS,rs)./rs);
    phiP     = rand(N,1)*2*pi;
    vxp      = vcircr.*cos(phiP);
    vyp      = vcircr.*sin(phiP);
    vx       = costheta.*(cos(phi).*vxp+sin(phi).*vyp);
    vy       = (-sin(phi).*vxp+cos(phi).*vyp);
    vz       = -sqrt(1-costheta.^2).*(cos(phi).*vxp+sin(phi).*vyp);
end






%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% RUN ODE SOL %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
Ncurrent=N;
MergeList=[];
tMin = 0;
tMax = TimeStepMassLoss;
dt= dtFactor*1e-5*2*pi*Rchar/sigmaS(find(rS>=Rchar,1,'first'));

inits=[X Y Z vx vy vz]; % initial conditions for each mass -- all in 1 column, 6 rows for ode solvers

options = odeset('reltol',1e-13*toleranceFactor,'abstol',1e-12*toleranceFactor,'Events',@MergerEvent); % tolerance for ode solver

epsilonTime = dt;
stops=0;

TTTg  = [];
XXXg  = [];
mtabg = [];

for kk = 1:round(SimTime/TimeStepMassLoss)
    %time_span=tMin:dt:tMax; % fixed mesh
    time_span=[tMin tMax];   % adaptive mesh, let ode45 figure it
    
    [TTT,XXX] = ode45(@nbodyCSoft,time_span,inits,options); %apparently this form is not as precise

    onelen = ones(length(TTT),1);
    mtab = onelen.*m';


    realMergeTidal=false;
    while TTT(end)+epsilonTime<tMax
        stops = stops + 1;
        disp(strcat('---Stop #',num2str(stops),' with MergerList ',num2str(MergeList(:,:))));

        mergerNum=1;
        if length(MergeList)>0
            disp(strcat('Merger! min/max:',num2str(min(m(MergeList(1,mergerNum)),m(MergeList(2,mergerNum)))),'/',num2str(max(m(MergeList(1,mergerNum)),m(MergeList(2,mergerNum))))))
            if min(m(MergeList(1,mergerNum)),m(MergeList(2,mergerNum))) > 1e-5 %make sure not bogus merger with itself
                realMergeTidal=true;
                disp('Actual merger!')
            end
        [initsnew, XXXnew]=StickParticles(MergeList(1,mergerNum),MergeList(2,mergerNum), XXX(end,:));
        elseif maxrc+0.01>=rTidal
            [initsnew, XXXnew]=CutTidalExit(XXX(end,:));
            disp('Exited tidal radius!')      
        end
        disp(strcat('Time at stop: ',num2str(TTT(end))));
        XXX(end,:) = XXXnew;

        TTT0 = TTT;
        XXX0 = XXX;
        tMin = TTT(end);
        time_span=tMin:dt:tMax;
        time_span = time_span(2:end);

        MergeList = []; 
        
        if length(time_span)==1
            time_span = [time_span (time_span+1e-7)]; % small bandaid if sims stop immediately, rare case
            disp('woops');
        end
        time_span = [time_span(1) time_span(end)]; % override with adaptive mesh
        
        [TTT,XXX] = ode45(@nbodyCSoft,time_span,initsnew,options);
        onelen = ones(length(TTT),1);
        mtab = [mtab' (onelen.*m')']';

        XXX = [XXX0' XXX']';
        TTT = [TTT0' TTT']';

        if realMergeTidal
            Ncurrent = Ncurrent - 1;
            realMergeTidal=false;
        end

    end
    
    m = mtab(end,:)' - MassLossPerStep;
    
    if HaloParameters.type == 1 % Stars-case: add mass to stars, modify properties
        MLRatio  = (massS(end)+sum(MassLossPerStep))/massS(end);
        massS    = MLRatio*massS;
        densityS = MLRatio*densityS;
        sigmaS   = MLRatio*sigmaS;
        phiS     = MLRatio*phiS;    
    end
    
    X=XXX(end,1:end/6);
    Y=XXX(end,(end/6+1):2*end/6);
    Z=XXX(end,(2*end/6+1):3*end/6);
        
    rt = sqrt(X.^2+Y.^2+Z.^2);

    [rt,sortIdx] = sort(rt,'ascend');
    % sort B using the sorting index
    m2 = m(sortIdx);
    mGCEnc = cumsum(m2);
    HaloMassInterp = interp1(rS,massS,rt);
    %mGCEncInterp = interp1(rt,mGCEnc,rS);
    %mGCEncInterp(isnan(mGCEncInterp))=0;
    idr=find(HaloMassInterp'-MGCMhaloRestrictFraction*mGCEnc>0,1,'first');
    %idr=find(massS-mGCEncInterp/2>0,1,'first');
    %size(massS-mGCEncInterp/2)
    rStallCandidate = rt(idr);
    %rStallCandidate = rS(idr)
    if ~isempty(idr) && rStallCandidate > HaloParameters.rStall
        rStall=rStallCandidate;
        %disp(strcat('changed stalling. frac: ',num2str(rStall/HaloParameters.rStall)))
    else
        rStall=HaloParameters.rStall;
        %disp(strcat('retained original stalling rad: '))
    end
    TauFudgeFac=((rS<rStall)*TauInStall+(rS>=rStall)*0)*TauFudgeFactorOverall;
        
    tMin = tMax+dt;
    tMax = tMin + TimeStepMassLoss;
    inits = [XXX(end,1:N)' XXX(end,(N+1):(2*N))' XXX(end,(2*N+1):(3*N))' XXX(end,(3*N+1):(4*N))' XXX(end,(4*N+1):(5*N))' XXX(end,(5*N+1):(6*N))'] ;
    XXXg = [XXXg' XXX']';
    TTTg = [TTTg' TTT']';
    mtabg = [mtabg' mtab']';
    
end


disp(strcat('total NaNs:',num2str(sum(sum(isnan(XXXg)))))) % Consistency. Should be zero. Typical reason for NaNs: GCs fly outside background function applicability (halo density, etc.)


SimStruct.mtab     = mtabg;
SimStruct.TTT      = TTTg;
SimStruct.XXX      = XXXg;
SimStruct.N        = N;
SimStruct.Ncurrent = Ncurrent;
SimStruct.Gn       = Gn;

end

%%
%% functions %%
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
function plotProjectedCDF(XX,YY,ZZ,saveFlag)
global m
minm=1e-8;
eligMass = m>minm;
% initial
rperp0A=(XX(1,eligMass).^2+YY(1,eligMass).^2).^0.5;
rperp0B=(ZZ(1,eligMass).^2+YY(1,eligMass).^2).^0.5;
rperp0C=(XX(1,eligMass).^2+ZZ(1,eligMass).^2).^0.5;
% final
rperp1A=(XX(end,eligMass).^2+YY(end,eligMass).^2).^0.5;
rperp1B=(ZZ(end,eligMass).^2+YY(end,eligMass).^2).^0.5;
rperp1C=(XX(end,eligMass).^2+ZZ(end,eligMass).^2).^0.5;

% group into empirical cdfs
[f0A,x0A]=ecdf(rperp0A);%,'Bounds','on')
[f0B,x0B]=ecdf(rperp0B);%,'Bounds','on')
[f0C,x0C]=ecdf(rperp0C);%,'Bounds','on')
[f1A,x1A]=ecdf(rperp1A);%,'Bounds','on')
[f1B,x1B]=ecdf(rperp1B);%,'Bounds','on')
[f1C,x1C]=ecdf(rperp1C);%,'Bounds','on')

figure
plot(x0A,f0A,'-','LineWidth',2)
hold on;
plot(x0B,f0B,'-','LineWidth',2)
plot(x0C,f0C,'-','LineWidth',2)
plot(x1A,f1A,'--','LineWidth',2)
plot(x1B,f1B,'--','LineWidth',2)
plot(x1C,f1C,'--','LineWidth',2)
grid on;
set(gca,'FontSize',15);
xlabel('$r_\perp$ [kpc]','interpreter','latex');
ylabel('CDF($r_\perp$)','interpreter','latex');
axis([0 4.5 0 1.05])
legend('Initial 1','Initial 2','Initial 3','Final 1','Final 2','Final 3', 'Location','NorthWest')      
if saveFlag
    saveas(gcf,'Plots/CDFperp.pdf')
end
end

function plotGCsInDiffMassBins(MassBins,XX,YY,ZZ,saveFlag)
global m
% get real data
binNum=length(MassBins)-1;
BinDat = ReadDataIntoMassBins(MassBins);

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
for ii=1:(length(MassBins)-1)
    %disp(strcat(num2str(MassBins(end-ii+1)),'and',num2str(MassBins(end-ii))))
    ind=find(msub>MassBins(end-ii));
    if ~isempty(ind)
        BinSim{(length(MassBins))-ii}=[rM1(ind) rM2(ind) rM3(ind)];
        msub(ind) = [];
        rM1(ind)   = [];
        rM2(ind)   = [];
        rM3(ind)   = [];    
    else
        BinSim{(length(MassBins))-ii}=[1000];
    end
end

co(1,:)=[0, 0.4470, 0.7410]; co(2,:)=[0.8500, 0.3250, 0.0980]; co(3,:)=[0.9290, 0.6940, 0.1250]; co(4,:)=[0.4940, 0.1840, 0.5560];
co(5,:)=[0.4660, 0.6740, 0.1880]; co(6,:)=[0.3010, 0.7450, 0.9330]; co(7,:)=[0.6350, 0.0780, 0.1840];
figure
hold on
%for jj=1:length(cdfBins)
%    [fc,xc]=ecdf(cdfBins{jj});
%    plot(xc,fc,'linewidth',2,'Color',co(jj,:))
%end
for jj=1:length(BinDat)
    [fc,xc]=ecdf(BinSim{jj});
    plot(xc,fc,'linewidth',2,'Color',co(jj,:))
    [fcDat,xcDat]=ecdf(BinDat{jj}(:,1));
    plot(xcDat,fcDat,'--','linewidth',2,'Color',co(jj,:))
end

grid on;
set(gca,'FontSize',15);
xlabel('$r_\perp$ [kpc]','interpreter','latex');
ylabel('CDF($r_\perp$)','interpreter','latex');
axis([0 4.5 0 1.05])
legendP=[];
for jj=1:length(MassBins)-1
    legendP = [legendP strcat("Bin ",num2str(jj), " sim") strcat("Bin ",num2str(jj), " data")];
end
legend(legendP, 'Location','SouthEast')
if ~isempty(MoreMassiveGCs)
    %disp('entered')
    title(strcat("Bins;",num2str(MassBins), " # larger GCs:",num2str(length(MoreMassiveGCs)),", max mass ",num2str(max(MoreMassiveGCs))),'interpreter','latex');
else
    title(strcat('Bins: ',num2str(MassBins)),'interpreter','latex');%; number of encounters ~ ', num2str(EncounterCount)))
end
if saveFlag
    saveas(gcf,'Plots/CDFperpBinned.pdf')
end
end

function plotCloseBound(TTT,tMax,XXX,critRadius,StepDiff)
numCloseArr=zeros(length(TTT),1);
numBoundArr=zeros(length(TTT),1);
numCloseAndBound=zeros(length(TTT),1);
for jj=1:StepDiff:length(TTT)
    coordsjj=XXX(jj,:)'; % use current state as initial conditions for next portion of orbit
    [CloseMatrix,NumClose,NumBound]=whosClose(coordsjj,critRadius);
    numCloseArr(jj) = NumClose;
    numBoundArr(jj) = NumBound;
    numCloseAndBound(jj) = HowManyCloseAndBound(coordsjj,critRadius);
end
figure
plot(TTT,numCloseArr)
grid on;
set(gca,'FontSize',15);
xlabel('t [Gyr]','interpreter','latex');
ylabel('GC pairs closeby','interpreter','latex');
axis([0 tMax 0 max(1,max(numCloseArr))+0.2])
%legend('Initial','Final', 'Location','NorthWest')
title(strcat('"close-by" radius = ',num2str(critRadius),' kpc'));%; number of encounters ~ ', num2str(EncounterCount)))

figure
plot(TTT,numBoundArr)
grid on;
set(gca,'FontSize',15);
xlabel('t [Gyr]','interpreter','latex');
ylabel('GC pairs bound, $E_{\rm rel} = \frac{\mu v_{\rm rel}^2}{2}+\Phi(r_{\rm rel}) < 0$','interpreter','latex');
axis([0 tMax 0 max(1,max(numBoundArr))+0.2])

figure
plot(TTT,numCloseAndBound)
grid on;
set(gca,'FontSize',15);
xlabel('t [Gyr]','interpreter','latex');
ylabel('Number of bound and close pairs','interpreter','latex');
axis([0 tMax 0 max(1,max(numCloseAndBound))+0.2])
end

function plotEnergy(TTT,CoarseTime,XX,YY,ZZ,VX,VY,VZ,saveFlag)
steps = length(TTT);
iDeltaT=find(TTT>=CoarseTime,1,'first');
timeCoarse = zeros(round(steps/iDeltaT),1);
energyCoarse = zeros(round(steps/iDeltaT),1);
c=0;
for i=1:iDeltaT:steps
    c=c+1;
    rTemp=[XX(i,:)' YY(i,:)' ZZ(i,:)'];
    timeCoarse(c) = TTT(i);
    energyCoarse(c) = getEnergy(VX(i,:),VY(i,:),VZ(i,:),rTemp);
end
figure
whitebg('white')
plot(timeCoarse,energyCoarse/abs(energyCoarse(1)),'linewidth',2)
grid on;
set(gca,'FontSize',15);
xlabel('t [Gyr]','interpreter','latex');
%ylabel('Energy [$10^5~M_\odot {\rm kpc}^2/{\rm Gyr}^2$]','interpreter','latex');
ylabel('-Energy/Initial Energy','interpreter','latex');
%axis([0 tMax 0 max(NumClosei)+0.5])
%legend('Initial','Final', 'Location','NorthWest')
%title(strcat('"close-by" radius = ',num2str(critRadius),' kpc'));%; number of encounters ~ ', num2str(EncounterCount)))

% display how much energy changed in console
r=[XX(1,:)' YY(1,:)' ZZ(1,:)'];
e0=getEnergy(VX(1,:),VY(1,:),VZ(1,:),r);
r=[XX(end,:)' YY(end,:)' ZZ(end,:)'];
e1=getEnergy(VX(end,:),VY(end,:),VZ(end,:),r);
disp(strcat('Initial/Final energy: ',num2str(round(e0/e1,3))))

if saveFlag
    saveas(gcf,'Plots/EnergyEvolution.pdf')
end
end

function plotMGCEvolution(TTT,mtab,tMax,N,Ncurrent,saveFlag)
figure
hold on;
for ii=1:N
    plot(TTT,mtab(:,ii),'linewidth',2)
end
grid on;
set(gca,'FontSize',15);
xlabel('t [Gyr]','interpreter','latex');
ylabel('M [$10^5~M_\odot$]','interpreter','latex');
axis([0 tMax 0 max(max(mtab))*1.1])
title(strcat('Initial GCs #: ',num2str(N),'. Final GC #: ',num2str(Ncurrent)))
if saveFlag
    saveas(gcf,'Plots/GCMassEvolution.pdf')
end

end

function [initsnew, XXXnew]=CutTidalExit(XXXlast)
global maxrc G Gn rS massS m rmerged realMergeTidal MassLossPerStep

N = length(XXXlast)/6;
i=1:N;
XX(i)=XXXlast(i);
YY(i)=XXXlast(i+N);
ZZ(i)=XXXlast(i+2*N);

VX(i)= XXXlast(i+3*N);
VY(i)= XXXlast(i+4*N);
VZ(i)= XXXlast(i+5*N);

rc = sqrt(XX.^2+YY.^2+ZZ.^2);
maxrcLocation = find(rc>=maxrc);

if m(maxrcLocation) >1e-5
    realMergeTidal=true;
    disp('Actual tidal disruption!')
end
m(maxrcLocation) = 1e-10;
MassLossPerStep(maxrcLocation) = 0;
%rmerged=10;
rmergedRand = rmerged*(0.8 + 0*0.33*rand(1,1));
vcircr   = sqrt(G*interp1(rS,massS,rmergedRand)./rmergedRand+Gn*sum(m)/rmergedRand);
phiP     = rand(1)*2*pi;
vxp      = vcircr.*cos(phiP);
vyp      = vcircr.*sin(phiP);
phi      = rand(1,1)*2*pi;
costheta = rand(1,1)*2-1;

XX(maxrcLocation) = rmergedRand*cos(phi).*sqrt(1-costheta.^2);
YY(maxrcLocation) = rmergedRand.*sin(phi).*sqrt(1-costheta.^2);
ZZ(maxrcLocation) = rmergedRand.*costheta;
VX(maxrcLocation) = costheta.*cos(phi).*vxp-sin(phi).*vyp;
VY(maxrcLocation) = sin(phi).*costheta.*vxp+cos(phi).*vyp;
VZ(maxrcLocation) = -sqrt(1-costheta.^2).*vxp;

initsnew = [XX' YY' ZZ' VX' VY' VZ'];
XXXnew   = [XX YY ZZ VX VY VZ];

end

function [newinits, XXXnew] = StickParticles(ii,jj, XXXlast)
global m rS massS G rmerged Gn MassLossPerStep
if m(ii)<m(jj)
    swap=ii;
    ii=jj;
    jj=swap;
end

N = length(XXXlast)/6;
i=1:N;
XX(i)=XXXlast(i);
YY(i)=XXXlast(i+N);
ZZ(i)=XXXlast(i+2*N);

VX(i)= XXXlast(i+3*N);
VY(i)= XXXlast(i+4*N);
VZ(i)= XXXlast(i+5*N);

XX(ii) = (m(ii)*XX(ii)+m(jj)*XX(jj))/(m(ii) + m(jj));
YY(ii) = (m(ii)*YY(ii)+m(jj)*YY(jj))/(m(ii) + m(jj));
ZZ(ii) = (m(ii)*ZZ(ii)+m(jj)*ZZ(jj))/(m(ii) + m(jj));

VX(ii) = (m(ii)*VX(ii)+m(jj)*VX(jj))/(m(ii) + m(jj));
VY(ii) = (m(ii)*VY(ii)+m(jj)*VY(jj))/(m(ii) + m(jj));
VZ(ii) = (m(ii)*VZ(ii)+m(jj)*VZ(jj))/(m(ii) + m(jj));

m(ii) = m(ii) + m(jj);
MassLossPerStep(ii) = MassLossPerStep(ii) + MassLossPerStep(jj);

m(jj) = 1e-10;
MassLossPerStep(jj) = 0;
%rmerged=10;
rmergedRand = rmerged*(0.8 + 0*0.33*rand(1,1));

vcircr   = sqrt(G*interp1(rS,massS,rmergedRand)./rmergedRand+Gn*sum(m)/rmergedRand);
phiP     = rand(1)*2*pi;
vxp      = vcircr.*cos(phiP);
vyp      = vcircr.*sin(phiP);
phi      = rand(1,1)*2*pi;
costheta = rand(1,1)*2-1;

XX(jj) = rmergedRand*cos(phi).*sqrt(1-costheta.^2);
YY(jj) = rmergedRand.*sin(phi).*sqrt(1-costheta.^2);
ZZ(jj) = rmergedRand.*costheta;
VX(jj) = costheta.*cos(phi).*vxp-sin(phi).*vyp;
VY(jj) = sin(phi).*costheta.*vxp+cos(phi).*vyp;
VZ(jj) = -sqrt(1-costheta.^2).*vxp;

newinits = [XX' YY' ZZ' VX' VY' VZ'];
XXXnew   = [XX YY ZZ VX VY VZ];
end

function [value, isterminal, direction] = MergerEvent(TTT, XXX)
    global m Gn MergeList critRadius maxrc rTidal epsilon
    %X=inits(:,1);
    %Y=inits(:,2);
    %Z=inits(:,3);
    %VX=inits(:,4);
    %VY=inits(:,5);
    %VZ=inits(:,6);
    N=length(XXX)/6;
    X=XXX(1:N);
    Y=XXX(N+1:2*N);
    Z=XXX(2*N+1:3*N);
    VX=XXX(3*N+1:4*N);
    VY=XXX(4*N+1:5*N);
    VZ=XXX(5*N+1:6*N);
    
    maxrc = max(sqrt(X.^2+Y.^2+Z.^2));
    
    DX=distance(X); % external function below
    DY=distance(Y);
    DZ=distance(Z);

    r=(DX.^2+DY.^2+DZ.^2).^0.5; % matrix of distances between all particles
    boolr=r<critRadius;
    
    %boolr = boolr - diag(diag(boolr)) % in case we want to set diagonal 0
    for ii=1:length(X)
        for jj=ii+1:length(X)
            if boolr(ii,jj)==1
                %E12 = 0.5*(m(ii)*m(jj)/(m(ii)+m(jj)))*((VX(ii)-VX(jj))^2+(VY(ii)-VY(jj))^2+(VZ(ii)-VZ(jj))^2) - Gn*m(ii)*m(jj)/r(ii,jj);
                E12 = 0.5*(m(ii)*m(jj)/(m(ii)+m(jj)))*((VX(ii)-VX(jj))^2+(VY(ii)-VY(jj))^2+(VZ(ii)-VZ(jj))^2) - Gn*m(ii)*m(jj)/(r(ii,jj)^2.11+(1.7*epsilon)^2.11)^(1/2.11);
                if E12<0
                    MergeList = [MergeList [ii, jj]'];
                end
            end
        end
    end
    %MergeList
    value      = [length(MergeList); (maxrc>rTidal)];%(coords(3) == 0.1);
    isterminal = [1; 1];   % Stop the integration
    direction  = [0; 0];
    %value       = [length(MergeList); maxrc
end

function en=getEnergy(VX,VY,VZ,r)
global m Gn G rS phiS epsilon
N=length(m);
V0=(VX.^2+VY.^2+VZ.^2).^0.5;
KE=sum(1/2*m'.*V0.^2,2); % total Kintetic Energy

% Potential Energy
PE=zeros(N,1);
for ii=1:N
    for jj=1:N
        if ii~=jj
            PE(ii)=PE(ii)-Gn*m(ii)*m(jj)./((r(jj,1)-r(ii,1))^2.11+(r(jj,2)-r(ii,2))^2.11+(r(jj,3)-r(ii,3))^2.11+(1.7*epsilon)^2.11)^(1/2.11)/2; % P.E. of mass i due to all j's
            %PE(ii)=PE(ii)-Gn*m(ii)*m(jj)./((r(jj,1)-r(ii,1))^2+(r(jj,2)-r(ii,2))^2+(r(jj,3)-r(ii,3))^2+epsilon^2)^0.5/2; % P.E. of mass i due to all j's
        end
    end
    PE(ii) = PE(ii) + m(ii)*interp1(rS,phiS,norm(r(ii,:)));%sqrt(r(ii,1)^2+r(ii,2)^2+r(ii,3)^2));
end
PE_tot=sum(PE);

en=KE+PE_tot;
end

function [boolr,num,numBound]=whosClose(inits,critRadius)
global m Gn
N=length(inits)/6;
X=inits(1:N);
Y=inits(N+1:2*N);
Z=inits(2*N+1:3*N);
VX=inits(3*N+1:4*N);
VY=inits(4*N+1:5*N);
VZ=inits(5*N+1:6*N);

DX=distance(X); % external function below
DY=distance(Y);
DZ=distance(Z);

r=(DX.^2+DY.^2+DZ.^2).^0.5; % matrix of distances between all particles

boolr=r<critRadius;

for jj=1:N
    boolr(jj,jj')=0;
end
num=sum(sum(boolr))/2;

rSq=distance(X).^2+distance(Y).^2+distance(Z).^2;
VSq=distance(VX).^2+distance(VY).^2+distance(VZ).^2;
muMat=m.*m'./(m+m');

ERel=0.5*muMat.*VSq - Gn*m.*m'./sqrt(rSq);
ERel(1:1+size(ERel,1):end) = 0;
numBound = sum(sum(ERel<0))/2;

end

function [numCloseAndBound]=HowManyCloseAndBound(inits,critRadius)
global m Gn epsilon
N=length(inits)/6;
X=inits(1:N);
Y=inits(N+1:2*N);
Z=inits(2*N+1:3*N);
VX=inits(3*N+1:4*N);
VY=inits(4*N+1:5*N);
VZ=inits(5*N+1:6*N);

DX=distance(X); % external function below
DY=distance(Y);
DZ=distance(Z);


r=(DX.^2+DY.^2+DZ.^2).^0.5; % matrix of distances between all particles
boolr=r<critRadius;

%boolr = boolr - diag(diag(boolr)) % in case we want to set diagonal 0
count = 0;
for ii=1:length(X)
    for jj=ii+1:length(X)
        if boolr(ii,jj)==1
            %E12 = 0.5*(m(ii)*m(jj)/(m(ii)+m(jj)))*((VX(ii)-VX(jj))^2+(VY(ii)-VY(jj))^2+(VZ(ii)-VZ(jj))^2) - Gn*m(ii)*m(jj)/r(ii,jj);
            E12 = 0.5*(m(ii)*m(jj)/(m(ii)+m(jj)))*((VX(ii)-VX(jj))^2+(VY(ii)-VY(jj))^2+(VZ(ii)-VZ(jj))^2) - Gn*m(ii)*m(jj)/(r(ii,jj)^2.11+(1.7*epsilon)^2.11)^(1/2.11);
            if E12<0
                count = count + 1;
            end
        end
    end
end

numCloseAndBound = count;
end

function [dist]=distance(X)
% this function calculates distance between two coordinates in a given vector X
dist = X-X.'; 
end
