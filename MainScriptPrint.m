%%
%MassBins = [exp(-0.2:0.8:4.5) 200]%[0.8 2 3 4.5 6 10 25 40 300];
clear

addpath(genpath('src'))

load(strcat('SimsBarMLS1.mat')) % load just to get some hyper parameters of sims. Be careful with this (in repurposing code)
%%
MassBins = [0.7 3 6 12 100];
cutRadius=3.8;
[rS0,TausS0,rB,TausB,rC,TausC] = getHaloRadialTau();

axisCommon = [0 20 0.8e-1 5];

Mx=0:1:20;

rTauEval=1.9;
TauBarEval = TausS0(find(rS0>rTauEval,1,'first'));
TauISOEval = TausB(find(rB>rTauEval,1,'first'));
TauNFWEval = TausC(find(rC>rTauEval,1,'first'));

rNFWCrit = rC(find(TausC>10,1,'first'));

[AvRProjBar,~] = getSersicMoments(nSersicGCsBar,GCSersicRadiusBar);
[AvRProjISO,~] = getSersicMoments(nSersicGCsISO,GCSersicRadiusISO);
[AvRProjNFW,~] = getSersicMoments(nSersicGCsNFW,GCSersicRadiusNFW);
%[AvRProjDense,~] = getSersicMoments(nSersicGCsDense,GCSersicRadiusDense);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Backbone Sims %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
modelName = 'BarMLS';modelNum  = '1'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['Stars. GCs obs. mass, $R_e^{(GC)}=' num2str(GCSersicRadiusBar) '$ kpc'],nSersicGCsBar,GCSersicRadiusBar,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjBar*exp(-(10/2/TauBarEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)
axis([0 20 0.8e-1 5])

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');

%%
modelName = 'ISOMLS';modelNum  = '1'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['Burkert. GCs of obs. mass, $R_e^{(GC)}=' num2str(GCSersicRadiusISO) '$ kpc'],nSersicGCsISO,GCSersicRadiusISO,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjISO*exp(-(10/2/TauISOEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');
%%
modelName = 'NFWMLS';modelNum  = '1'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['NFW. GCs of obs. mass, $R_e^{(GC)}=' num2str(GCSersicRadiusNFW) '$ kpc'],nSersicGCsNFW,GCSersicRadiusNFW,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjNFW*exp(-(10/2/TauNFWEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Backbone Sim: Almost obs %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
modelName = 'NFWMLSAO';modelNum  = '1'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['NFW. GCs almost as obs., $R_e^{(GC)}=' num2str(GCSersicRadiusNFW) '$ kpc'],nSersicGCsNFW,GCSersicRadiusNFW,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjNFW*exp(-(10/2/TauNFWEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)
axis([0 26 0.2e-1 5])

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Extra sim: denser GC sim %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
MassBins = [0.7 3 6 12 30 100];

modelName = 'NFWMLSDG';modelNum  = '1'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['NFW. GCs as obs., $R_e^{(GC)}=' num2str(GCSersicRadiusNFW) '$ kpc'],nSersicGCsNFW,GCSersicRadiusNFW,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjNFW*exp(-(10/2/TauNFWEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)
axis([0 45 0.06e-1 5])

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');

%%
modelName = 'NFWMLSDG';modelNum  = '1UR'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['NFW. GCs as obs., $R_e^{(GC)}=' num2str(GCSersicRadiusNFW) '$ kpc, DF U.R.'],nSersicGCsNFW,GCSersicRadiusNFW,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjNFW*exp(-(10/2/TauNFWEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)
axis([0 45 0.06e-1 5])

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');

%%
modelName = 'ISOMLSDG';modelNum  = '1'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['ISO. GCs as obs., $R_e^{(GC)}=' num2str(GCSersicRadiusISO) '$ kpc'],nSersicGCsISO,GCSersicRadiusISO,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjISO*exp(-(10/2/TauISOEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)
%axis([0 45 0.06e-1 5])

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');
%%
modelName = 'ISOMLSDG';modelNum  = '1UR'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['ISO. GCs as obs., $R_e^{(GC)}=' num2str(GCSersicRadiusISO) '$ kpc, DF U.R.'],nSersicGCsISO,GCSersicRadiusISO,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjISO*exp(-(10/2/TauISOEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)
%axis([0 45 0.06e-1 5])

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Unrestricted Sims %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
modelName = 'BarMLSUR';modelNum  = '1'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['Stars. GCs obs. mass, $R_e^{(GC)}=' num2str(GCSersicRadiusBar) '$ kpc, DF U.R.'],nSersicGCsBar,GCSersicRadiusBar,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjBar*exp(-(10/2/TauBarEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)
axis([0 27 0.8e-1 5])

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');

%%
modelName = 'ISOMLSUR';modelNum  = '1'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['Burkert. GCs of obs. mass, $R_e^{(GC)}=' num2str(GCSersicRadiusISO) '$ kpc, DF U.R.'],nSersicGCsISO,GCSersicRadiusISO,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjISO*exp(-(10/2/TauISOEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');
%%
modelName = 'NFWMLSUR';modelNum  = '1'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['NFW. GCs of obs. mass, $R_e^{(GC)}=' num2str(GCSersicRadiusNFW) '$ kpc, DF U.R.'],nSersicGCsNFW,GCSersicRadiusNFW,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjNFW*exp(-(10/2/TauNFWEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Extra sim: Longer GC sim (UR) %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NFW
MassBins = [0.7 3 6 12 40 60 100];

modelName = 'NFWMLSLongUR';modelNum  = '1'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))

getLumRadVsMassSims3TwoCol(MassBins,St1A,['NFW. GCs as obs., $R_e^{(GC)}=' num2str(GCSersicRadiusNFW) '$ kpc, 20 Gyr'],nSersicGCsNFW,GCSersicRadiusNFW,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjNFW*exp(-(20/2/TauNFWEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)
axis([0 70 0.06e-1 5])

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');

%% ISO
MassBins = [0.7 3 6 12 40 60 100];

modelName = 'ISOMLSLongUR';modelNum  = '1'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))

getLumRadVsMassSims3TwoCol(MassBins,St1A,['ISO. GCs as obs., $R_e^{(GC)}=' num2str(GCSersicRadiusISO) '$ kpc, 20 Gyr'],nSersicGCsISO,GCSersicRadiusISO,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjISO*exp(-(20/2/TauISOEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)
axis([0 45 0.06e-1 5])

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Extra sim: Longer GC sim - less merger %%%%SimsNFWMLSLongTrem1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NFW
MassBins = [0.7 3 6 12 40 60 100];

modelName = 'NFWMLSLongTrem';modelNum  = '1'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))

getLumRadVsMassSims3TwoCol(MassBins,St1A,['NFW. GCs as obs., $R_e^{(GC)}=' num2str(GCSersicRadiusNFW) '$ kpc, 20 Gyr'],nSersicGCsNFW,GCSersicRadiusNFW,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjNFW*exp(-(20/2/TauNFWEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)
axis([0 70 0.06e-1 5])

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');

%% ISO
MassBins = [0.7 3 6 12 40 60 100];

modelName = 'ISOMLSLongTrem';modelNum  = '1'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))

getLumRadVsMassSims3TwoCol(MassBins,St1A,['ISO. GCs as obs., $R_e^{(GC)}=' num2str(GCSersicRadiusISO) '$ kpc, 20 Gyr'],nSersicGCsISO,GCSersicRadiusISO,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjISO*exp(-(20/2/TauISOEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)
axis([0 45 0.06e-1 5])

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Low Eps %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
modelName = 'BarMLS';modelNum  = '1EpsLow'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['Stars. GCs obs. mass, $R_e^{(GC)}=' num2str(GCSersicRadiusBar) '$ kpc, $\epsilon=4~$pc'],nSersicGCsBar,GCSersicRadiusBar,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjBar*exp(-(10/2/TauBarEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)
%axis([0 20 0.8e-1 5])

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');

%%
modelName = 'ISOMLS';modelNum  = '1EpsLow'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['Burkert. GCs of obs. mass, $R_e^{(GC)}=' num2str(GCSersicRadiusISO) '$ kpc, $\epsilon=4~$pc'],nSersicGCsISO,GCSersicRadiusISO,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjISO*exp(-(10/2/TauISOEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');
%%
modelName = 'NFWMLS';modelNum  = '1EpsLow'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['NFW. GCs of obs. mass, $R_e^{(GC)}=' num2str(GCSersicRadiusNFW) '$ kpc, $\epsilon=4~$pc'],nSersicGCsNFW,GCSersicRadiusNFW,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjNFW*exp(-(10/2/TauNFWEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% High Eps %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
modelName = 'BarMLS';modelNum  = '1EpsHigh'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['Stars. GCs obs. mass, $R_e^{(GC)}=' num2str(GCSersicRadiusBar) '$ kpc, $\epsilon=12~$pc'],nSersicGCsBar,GCSersicRadiusBar,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjBar*exp(-(10/2/TauBarEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)
%axis([0 20 0.8e-1 5])

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');

%%
modelName = 'ISOMLS';modelNum  = '1EpsHigh'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['Burkert. GCs of obs. mass, $R_e^{(GC)}=' num2str(GCSersicRadiusISO) '$ kpc, $\epsilon=12~$pc'],nSersicGCsISO,GCSersicRadiusISO,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjISO*exp(-(10/2/TauISOEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');
%%
modelName = 'NFWMLS';modelNum  = '1EpsHigh'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['NFW. GCs of obs. mass, $R_e^{(GC)}=' num2str(GCSersicRadiusNFW) '$ kpc, $\epsilon=12~$pc'],nSersicGCsNFW,GCSersicRadiusNFW,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjNFW*exp(-(10/2/TauNFWEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Low Merger radius %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
modelName = 'BarMLS';modelNum  = '1MrgLow'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['Stars. GCs obs. mass, $R_e^{(GC)}=' num2str(GCSersicRadiusBar) '$ kpc, $r_{\rm merger}=10~$pc'],nSersicGCsBar,GCSersicRadiusBar,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjBar*exp(-(10/2/TauBarEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)
%axis([0 20 0.8e-1 5])

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');

%%
modelName = 'ISOMLS';modelNum  = '1MrgLow'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['Burkert. GCs of obs. mass, $R_e^{(GC)}=' num2str(GCSersicRadiusISO) '$ kpc, $r_{\rm merger}=10~$pc'],nSersicGCsISO,GCSersicRadiusISO,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjISO*exp(-(10/2/TauISOEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');
%%
modelName = 'NFWMLS';modelNum  = '1MrgLow'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['NFW. GCs of obs. mass, $R_e^{(GC)}=' num2str(GCSersicRadiusNFW) '$ kpc, $r_{\rm merger}=10~$pc'],nSersicGCsNFW,GCSersicRadiusNFW,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjNFW*exp(-(10/2/TauNFWEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% High Merger radius %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
modelName = 'BarMLS';modelNum  = '1MrgHigh'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['Stars. GCs obs. mass, $R_e^{(GC)}=' num2str(GCSersicRadiusBar) '$ kpc, $r_{\rm merger}=35~$pc'],nSersicGCsBar,GCSersicRadiusBar,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjBar*exp(-(10/2/TauBarEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)
axis([0 24 0.8e-1 5])

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');

%%
modelName = 'ISOMLS';modelNum  = '1MrgHigh'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['Burkert. GCs of obs. mass, $R_e^{(GC)}=' num2str(GCSersicRadiusISO) '$ kpc, $r_{\rm merger}=35~$pc'],nSersicGCsISO,GCSersicRadiusISO,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjISO*exp(-(10/2/TauISOEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');
%%
modelName = 'NFWMLS';modelNum  = '1MrgHigh'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TwoCol(MassBins,St1A,['NFW. GCs of obs. mass, $R_e^{(GC)}=' num2str(GCSersicRadiusNFW) '$ kpc, $r_{\rm merger}=35~$pc'],nSersicGCsNFW,GCSersicRadiusNFW,cutRadius,0,saveName,true);

subplot(2,1,1)
plot(Mx,AvRProjNFW*exp(-(10/2/TauNFWEval)*(Mx/5).^1),'--','linewidth',1.3)
legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Backbone Sims Only top panel! %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
modelName = 'BarMLS';modelNum  = '1'; saveName=strcat(modelName,modelNum,'s1ATop');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TopPanel(MassBins,St1A,['\bf Stars. $R_e^{(GC)}=' num2str(GCSersicRadiusBar) '$ kpc'],nSersicGCsBar,GCSersicRadiusBar,cutRadius,0,saveName,true);

%subplot(2,1,1)
%plot(Mx,AvRProjBar*exp(-(10/2/TauBarEval)*(Mx/5).^1),'--','linewidth',1.3)
%legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
legend('Data','Simulation','Initial condition','Poor theory control','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)
axis([0 20 0.8e-1 5])

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');

%%
modelName = 'ISOMLS';modelNum  = '1'; saveName=strcat(modelName,modelNum,'s1ATop');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TopPanel(MassBins,St1A,['\bf Burkert. $R_e^{(GC)}=' num2str(GCSersicRadiusISO) '$ kpc'],nSersicGCsISO,GCSersicRadiusISO,cutRadius,0,saveName,true);

%subplot(2,1,1)
%plot(Mx,AvRProjISO*exp(-(10/2/TauISOEval)*(Mx/5).^1),'--','linewidth',1.3)
%legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
legend('Data','Simulation','Initial condition','Poor theory control','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');
%%
modelName = 'NFWMLS';modelNum  = '1'; saveName=strcat(modelName,modelNum,'s1ATop');
load(strcat('Sims',modelName,modelNum,'.mat'))
getLumRadVsMassSims3TopPanel(MassBins,St1A,['\bf NFW. $R_e^{(GC)}=' num2str(GCSersicRadiusNFW) '$ kpc'],nSersicGCsNFW,GCSersicRadiusNFW,cutRadius,0,saveName,true);

%subplot(2,1,1)
%plot(Mx,AvRProjNFW*exp(-(10/2/TauNFWEval)*(Mx/5).^1),'--','linewidth',1.3)
%legend('Data','Simulation','Initial condition','Poor theory control','Eq. (2) approx.','Location','SouthWest')%,'interpreter','latex'
legend('Data','Simulation','Initial condition','Poor theory control','Location','SouthWest')%,'interpreter','latex'
axis(axisCommon)

saveas(gcf,strcat('Plots/',saveName,'TP.eps'),'epsc');

