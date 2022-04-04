% This script runs the main simulations Bar, Danieli & Blum, arXiv:2202.10179
%
% It stars by defining the parameters for simulations in different halos,
% followed by a suite of ~40 runs for each halo configuration.
% 
% To generate the plots in the paper (Fig. 5, 6 and appendices), refer to
% MainScriptPrint.m
%
% To generate a video of a sample simulation, refer to MainScriptVideo.m

clear all;
addpath(genpath('src'))
addpath(genpath('data'))

%% set default models and parameters

% Define different halos via "config" structures
HaloParamBar.type = 1; % 1-Sersic
HaloParamBar.nSersic = 0.61;
HaloParamBar.MSersic = 1.2e3;
HaloParamBar.rSersic = 1.9;
HaloParamBar.rStall = (HaloParamBar.rSersic)*0.3;

HaloParamNFW.type = 2; % 2-NFW
HaloParamNFW.NFWRs = 6;
HaloParamNFW.NFWc = 6;
HaloParamNFW.rStall = 0.01;

HaloParamISO.type = 3; % 3 - Burkert
HaloParamISO.BurRs = 2;
HaloParamISO.Burc = 15;
HaloParamISO.rStall = (HaloParamISO.BurRs)*0.3; %*0.3

% GCs initial 
GCSersicRadiusBar = 4.5; nSersicGCsBar = 0.61; % baryon
GCSersicRadiusNFW = 2.5; nSersicGCsNFW = 0.61;
GCSersicRadiusISO = 2.6; nSersicGCsISO = 0.61; 

% some other options
GnFrac=1;
TauFudgeFactorOverll=1;
SimTime=10.05;

TidalRadius=20; %where to cut GCs
SofteningRadKpc=0.007; % Plummer-softening of GCs
IsoVsCircFlag=0;       % 0: isotropic dist, 1: circ dist

tolFactor=3e6; % in comparison to 1e-13/1e-12 rel/abs or so built-in the code

NUMExamples=1; % simulations number in each 

deltaTCoarseTime = 0.005; % time-step after saving 
dtFactor = 0.5; 

MergeRadius = 0.02;
MGCMhaloRestrictFraction = 1/2;
TotalLossFraction = (1/3)/10; % third mass loss over 10Gyr

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Paper backbone (Fig. 5) %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NFW halo 1
HaloParam = HaloParamNFW;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusNFW,nSersicGCsNFW,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsNFWMLS1.mat')
clear St1A
%% ISO halo 1
HaloParam = HaloParamISO;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusISO,nSersicGCsISO,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsISOMLS1.mat')
clear St1A
%% BAR halo 1
HaloParam = HaloParamBar;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusBar,nSersicGCsBar,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsBarMLS1.mat')
clear St1A


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Paper backbone (Fig. 6) Almost Obs. %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NFW halo ALMOST OBSERVED
HaloParam = HaloParamNFW;
[St1A]=runSuiteSameStarAlmostObs(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusNFW,nSersicGCsNFW,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsNFWMLSAO1.mat')
clear St1A

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% "Denser" GC distribution %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NFW halo Denser
GCSersicRadiusNFWDef = GCSersicRadiusNFW;
GCSersicRadiusNFW = 1.9;
HaloParam = HaloParamNFW;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusNFW,nSersicGCsNFW,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsNFWMLSDG1.mat') % DG = Dense GCs
clear St1A
GCSersicRadiusNFW = GCSersicRadiusNFWDef;
%% ISO halo Denser
GCSersicRadiusISODef = GCSersicRadiusISO;
GCSersicRadiusISO = 1.9; 
HaloParam = HaloParamISO;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusISO,nSersicGCsISO,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsISOMLSDG1.mat') % DG = Dense GCs
clear St1A
GCSersicRadiusISO = GCSersicRadiusISODef;

%% UR
HaloParamISODef = HaloParamISO;
HaloParamNFWDef = HaloParamNFW;
MGCMhaloRestrictFractionDef = MGCMhaloRestrictFraction;

HaloParamBar.rStall = 0.01;
HaloParamISO.rStall = 0.01; %*0.3
HaloParamNFW.rStall = 0.01;%sqrt(30/(2*pi*0.001359*(200/3)*(HaloParamNFW.NFWc)^3*(log(1+HaloParamNFW.NFWc)-HaloParamNFW.NFWc/(1+HaloParamNFW.NFWc))^(-1)*HaloParamNFW.NFWRs)); % roughly radius enclosing 30e5 Msol
MGCMhaloRestrictFraction = 0.01;

%% NFW halo Denser UR
GCSersicRadiusNFWDef = GCSersicRadiusNFW;
GCSersicRadiusNFW = 1.9;
HaloParam = HaloParamNFW;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusNFW,nSersicGCsNFW,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsNFWMLSDG1UR.mat') % DG = Dense GCs
clear St1A
GCSersicRadiusNFW = GCSersicRadiusNFWDef;
%% ISO halo Denser UR
GCSersicRadiusISODef = GCSersicRadiusISO;
GCSersicRadiusISO = 1.9; 
HaloParam = HaloParamISO;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusISO,nSersicGCsISO,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsISOMLSDG1UR.mat') % DG = Dense GCs
clear St1A
GCSersicRadiusISO = GCSersicRadiusISODef;
%%
HaloParamISO = HaloParamISODef;
HaloParamNFW = HaloParamNFWDef;
MGCMhaloRestrictFraction = MGCMhaloRestrictFractionDef;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% "Later" GC distribution %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SimTimeDef = SimTime;
SimTime = 20;
MGCMhaloRestrictFractionDef = MGCMhaloRestrictFraction;
MGCMhaloRestrictFraction = 0.01;
%% NFW halo 1
HaloParam = HaloParamNFW;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusNFW,nSersicGCsNFW,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsNFWMLSLongUR1.mat')
clear St1A
%% ISO halo 1
HaloParam = HaloParamISO;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusISO,nSersicGCsISO,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsISOMLSLongUR1.mat')
clear St1A
%%
MGCMhaloRestrictFraction = MGCMhaloRestrictFractionDef;
SimTime = SimTimeDef;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% "Later" GC distribution No DF No Mergers %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SimTimeDef = SimTime;
SimTime = 15;
MGCMhaloRestrictFractionDef = MGCMhaloRestrictFraction;
MGCMhaloRestrictFraction = 0.01;
MergeRadiusDef = MergeRadius;
MergeRadius = 0.02;
GnFracDef=GnFrac;
GnFrac = 0.001;
%TotalLossFractionDef = TotalLossFraction;
%TotalLossFraction = 1e-6;
%% NFW halo 1
HaloParam = HaloParamNFW;
[St1A]=runSuiteSameStar(deltaTCoarseTime,10,HaloParam,GCSersicRadiusNFW,nSersicGCsNFW,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsNFWMLSLongTrem1.mat')
clear St1A
%% ISO halo 1
HaloParam = HaloParamISO;
[St1A]=runSuiteSameStar(deltaTCoarseTime,3,HaloParam,GCSersicRadiusISO,nSersicGCsISO,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsISOMLSLongTrem1.mat')
clear St1A
%%
GnFrac=GnFracDef;
MGCMhaloRestrictFraction = MGCMhaloRestrictFractionDef;
SimTime = SimTimeDef;
MergeRadius = MergeRadiusDef;
%TotalLossFraction = TotalLossFractionDef;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Unrestricted DF to origin %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HaloParamBarDef = HaloParamBar;
HaloParamISODef = HaloParamISO;
HaloParamNFWDef = HaloParamNFW;
MGCMhaloRestrictFractionDef = MGCMhaloRestrictFraction;

HaloParamBar.rStall = 0.01;
HaloParamISO.rStall = 0.01; %*0.3
HaloParamNFW.rStall = 0.01;%sqrt(30/(2*pi*0.001359*(200/3)*(HaloParamNFW.NFWc)^3*(log(1+HaloParamNFW.NFWc)-HaloParamNFW.NFWc/(1+HaloParamNFW.NFWc))^(-1)*HaloParamNFW.NFWRs)); % roughly radius enclosing 30e5 Msol
MGCMhaloRestrictFraction = 0.01;

%% NFW halo 1
HaloParam = HaloParamNFW;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusNFW,nSersicGCsNFW,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsNFWMLSUR1.mat')
clear St1A
%% ISO halo 1
HaloParam = HaloParamISO;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusISO,nSersicGCsISO,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsISOMLSUR1.mat')
clear St1A
%% BAR halo 1
HaloParam = HaloParamBar;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusBar,nSersicGCsBar,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsBarMLSUR1.mat')
clear St1A
%%
HaloParamBar = HaloParamBarDef;
HaloParamISO = HaloParamISODef;
HaloParamNFW = HaloParamNFWDef;
MGCMhaloRestrictFraction = MGCMhaloRestrictFractionDef;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Appendix tests (Fig. 5) %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Low epsilon
SofteningRadKpcDef = SofteningRadKpc;
SofteningRadKpc = 0.004;
%% NFW halo 1
HaloParam = HaloParamNFW;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusNFW,nSersicGCsNFW,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsNFWMLS1EpsLow.mat')
clear St1A
%% ISO halo 1
HaloParam = HaloParamISO;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusISO,nSersicGCsISO,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsISOMLS1EpsLow.mat')
clear St1A
%% BAR halo 1
HaloParam = HaloParamBar;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusBar,nSersicGCsBar,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsBarMLS1EpsLow.mat')
clear St1A
SofteningRadKpc = SofteningRadKpcDef;
%% High epsilon
SofteningRadKpcDef = SofteningRadKpc;
SofteningRadKpc = 0.012;
%% NFW halo 1
HaloParam = HaloParamNFW;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusNFW,nSersicGCsNFW,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsNFWMLS1EpsHigh.mat')
clear St1A
%% ISO halo 1
HaloParam = HaloParamISO;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusISO,nSersicGCsISO,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsISOMLS1EpsHigh.mat')
clear St1A
%% BAR halo 1
HaloParam = HaloParamBar;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusBar,nSersicGCsBar,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsBarMLS1EpsHigh.mat')
clear St1A
SofteningRadKpc = SofteningRadKpcDef;


%% Low merger radius
MergeRadiusDef = MergeRadius;
MergeRadius = 0.01;
%% NFW halo 1
HaloParam = HaloParamNFW;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusNFW,nSersicGCsNFW,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsNFWMLS1MrgLow.mat')
clear St1A
%% ISO halo 1
HaloParam = HaloParamISO;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusISO,nSersicGCsISO,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsISOMLS1MrgLow.mat')
clear St1A
%% BAR halo 1
HaloParam = HaloParamBar;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusBar,nSersicGCsBar,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsBarMLS1MrgLow.mat')
clear St1A
MergeRadius = MergeRadiusDef;
%% High merger radius
MergeRadiusDef = MergeRadius;
MergeRadius = 0.035;
%% NFW halo 1
HaloParam = HaloParamNFW;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusNFW,nSersicGCsNFW,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsNFWMLS1MrgHigh.mat')
clear St1A
%% ISO halo 1
HaloParam = HaloParamISO;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusISO,nSersicGCsISO,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsISOMLS1MrgHigh.mat')
clear St1A
%% BAR halo 1
HaloParam = HaloParamBar;
[St1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadiusBar,nSersicGCsBar,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction);
save('SimsBarMLS1MrgHigh.mat')
clear St1A
MergeRadius = MergeRadiusDef;

%% FUNCTIONS

function [Struct1A]=runSuiteSameStar(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadius,nSersicGCs,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction)
% run with observed GC distribution. 

massScaleIn=ReadDataIntoMassSpec(); 
Nin = length(massScaleIn);
MassNormMean=-1; MassNormVar=0.3; % not relevant when massScaleIn is array

rng(0);
for ii=1:NUMExamples
    SimStructT = RunSimMassLossStall(HaloParam,GCSersicRadius,nSersicGCs,Nin,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,massScaleIn,MassNormMean,MassNormVar,MergeRadius,TidalRadius,tolFactor,dtFactor,MGCMhaloRestrictFraction,TotalLossFraction);
    Struct1A{ii} = ReduceStruct(SimStructT,deltaTCoarseTime,MassNormMean,MassNormVar);
    clear SimStructT
end

end

function [Struct1A]=runSuiteSameStarAlmostObs(deltaTCoarseTime,NUMExamples,HaloParam,GCSersicRadius,nSersicGCs,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,TidalRadius,tolFactor,dtFactor,MergeRadius,MGCMhaloRestrictFraction,TotalLossFraction)
% run with ALMOST observed GC distribution. GCIMFSt irrelevant

massScaleIn = getGCMassSpecAlmostObserved();
Nin = length(massScaleIn);
MassNormMean=-3; MassNormVar=0.3; % not relevant when massScaleIn is array except as flag
rng(0);
for ii=1:NUMExamples
    SimStructT = RunSimMassLossStall(HaloParam,GCSersicRadius,nSersicGCs,Nin,SimTime,GnFrac,TauFudgeFactorOverll,IsoVsCircFlag,SofteningRadKpc,massScaleIn,MassNormMean,MassNormVar,MergeRadius,TidalRadius,tolFactor,dtFactor,MGCMhaloRestrictFraction,TotalLossFraction);
    Struct1A{ii} = ReduceStruct(SimStructT,deltaTCoarseTime,MassNormMean,MassNormVar);
    clear SimStructT
end

end


function sx2 = ReduceStruct(sx,deltaT,M0,sigma)
%deltaT=0.01;
NT=round(sx.TTT(end)/deltaT);
TTTind=zeros(length(NT),1);
TTTind(1)=1;
for k=2:NT
    TTTind(k)=find(sx.TTT>=deltaT*(k-1),1,'first');
end
sx2 = sx;

sx2.TTT = sx.TTT(TTTind);
sx2.XXX = sx.XXX(TTTind,:);
sx2.mtab = sx.mtab(TTTind,:);

sx2.MeanInitial = M0;
sx2.SigmaInitial = sigma;
end


