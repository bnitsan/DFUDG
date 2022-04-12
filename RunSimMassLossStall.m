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
cutSersicId = find(rSGCs > rTidal,1,'first'); % rarely, for very diffuse distributions, the sampling may generate GCs very far. We account for that by cutting the CDF at the "tidal" radius
SersicCDF = massSGCs/massSGCs(cutSersicId);
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
        if rr(1) >= rc0(kk)
            ix = 1;
        elseif rr(end)<= rc0(kk)
            ix = length(rr);
        else
            ix=find(rr>rc0(kk),1,'first');
        end
        disp([num2str(rc0(kk)) ',' num2str(ix)])

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

% every loop the mass of GCs is reduced as per mass loss rate
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

    % here we enforce the condition that DF is ineffective for 
    % Mhalo-frac*MGCenclosed < 0.
    [rt,sortIdx] = sort(rt,'ascend');
    m2 = m(sortIdx);
    mGCEnc = cumsum(m2);
    HaloMassInterp = interp1(rS,massS,rt);
    idr=find(HaloMassInterp'-MGCMhaloRestrictFraction*mGCEnc>0,1,'first');
    rStallCandidate = rt(idr);

    if ~isempty(idr) && rStallCandidate > HaloParameters.rStall
        rStall=rStallCandidate;
    else
        rStall=HaloParameters.rStall;
    end
    TauFudgeFac=((rS<rStall)*TauInStall+(rS>=rStall)*0)*TauFudgeFactorOverall;
        
    % rebuild simulation structures from this loop
    tMin = tMax+dt;
    tMax = tMin + TimeStepMassLoss;
    inits = [XXX(end,1:N)' XXX(end,(N+1):(2*N))' XXX(end,(2*N+1):(3*N))' XXX(end,(3*N+1):(4*N))' XXX(end,(4*N+1):(5*N))' XXX(end,(5*N+1):(6*N))'] ;
    XXXg = [XXXg' XXX']';
    TTTg = [TTTg' TTT']';
    mtabg = [mtabg' mtab']';
    
end


disp(strcat('total NaNs:',num2str(sum(sum(isnan(XXXg)))))) % Consistency. Should be zero. Typical reason for NaNs: GCs fly outside background function applicability (halo density, etc.)

% Construct a struct that contains output of simulation
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

function [initsnew, XXXnew]=CutTidalExit(XXXlast)
% this function cuts a GC beyond the tidal radius. Gets XXX structure of
% the simulations, returns new values for the continuation of the
% simulation
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
% a function that implements mergers by "sticking" GCs ii and jj in XXXlast
% of the simulation. Returns new phase-space structures to continue the simulation
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
% a function that detects a merger event or tidal disruption events
    global m Gn MergeList critRadius maxrc rTidal epsilon

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
    
    for ii=1:length(X)
        for jj=ii+1:length(X)
            if boolr(ii,jj)==1
                E12 = 0.5*(m(ii)*m(jj)/(m(ii)+m(jj)))*((VX(ii)-VX(jj))^2+(VY(ii)-VY(jj))^2+(VZ(ii)-VZ(jj))^2) - Gn*m(ii)*m(jj)/(r(ii,jj)^2.11+(1.7*epsilon)^2.11)^(1/2.11);
                if E12<0
                    MergeList = [MergeList [ii, jj]'];
                end
            end
        end
    end
    value      = [length(MergeList); (maxrc>rTidal)];%(coords(3) == 0.1);
    isterminal = [1; 1];   % Stop the integration
    direction  = [0; 0];
end

function [dist]=distance(X)
% this function calculates distance between two coordinates in a given vector X
dist = X-X.'; 
end
