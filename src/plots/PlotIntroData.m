clear
%%

MassToLightRatio=1.6;

[rperp,MGC] = ReadDataIntoMassSpecAndRad();

mbins3=exp(log(0.7):0.8:log(30));

[AvgMass1,AvgRProj1,sizePoint1,NumBin1]=getMomentsMR([0.7 3 6 12 25]);
[AvgMass2,AvgRProj2,sizePoint2,NumBin2]=getMomentsMR(exp(log(0.7):0.9:log(30)));
[AvgMass3,AvgRProj3,sizePoint3,NumBin3]=getMomentsMR(exp(log(0.7):0.8:log(30)));

% in case last bin aligns in two binnings, shift a little
AvgRProj1(end) = AvgRProj1(end)+0.01;
AvgRProj2(end) = AvgRProj2(end)-0.01;

c2=[0.8500, 0.3250, 0.0980];
c1=[0.4940, 0.1840, 0.5560]	;
figure
scatter(MGC/MassToLightRatio,rperp,150)
hold on;
%scatter(AvgMass1,AvgRProj1,200,c1);%sizePoint1);
%scatter(AvgMass2,AvgRProj2,sizePoint2,'filled');
scatter(AvgMass3/MassToLightRatio,AvgRProj3,250,c2);%sizePoint3);

%text(AvgMass1-0.23,AvgRProj1+0.01,num2str(NumBin1),'Color',c1,'FontSize',14)
text(AvgMass3/MassToLightRatio-0.18,AvgRProj3+0.02,num2str(NumBin3),'Color',c2,'FontSize',14)
%scatter(AvgMass4,AvgRProj4,sizePoint4,'filled');
%set(gca,'YScale','log')
%axis([0.0 max(max(SimsAvgMass),max(DatAvgMass))*1.2 0 4])
xlabel('L bins [$10^5~L_\odot$]','interpreter','latex');
ylabel('$\left<r_\perp\right>$ [kpc]','interpreter','latex')

axis([0 10 0 4])
grid on;
set(gca,'FontSize',15);
%legend('All data','Binning 1','Binning 2')
%title(titleFig)

for ii=2:(length(mbins3)-1)
    xline(mbins3(ii)/MassToLightRatio,'--')
end
    
legend('All data','Luminosity-binned data')

%saveas(gcf,strcat('Plots/IntroCleanData.eps'),'epsc');

%%
figure
scatter(AvgMass1,AvgRProj1-0.005,sizePoint1,'filled');
hold on;
scatter(AvgMass2,AvgRProj2+0.005,sizePoint2,'filled');
scatter(AvgMass3,AvgRProj3,sizePoint3,'filled');
Mx=0:0.2:20;
plot(Mx,3.3*exp(-(10/2/5.0)*(Mx/5).^1),'linewidth',2)
set(gca,'YScale','log')

%axis([0.0 max(max(SimsAvgMass),max(DatAvgMass))*1.2 0 4])
xlabel('M bins [$10^5~M_\odot$]','interpreter','latex');
ylabel('$\left<r_\perp\right>$ [kpc]','interpreter','latex')
grid on;
set(gca,'FontSize',15);
%legend('Data','Simple dynamical friction model')
legend('Data Binning 1','Data Binning 2','Data Binning 3','Simple dynamical friction model')
%title(titleFig)
axis([0.5 16 0.1 4])
saveas(gcf,strcat('Plots/IntroCleanDataWithFit.eps'),'epsc');

%%
c1=[0, 0.4470, 0.7410];
c2=[0.8500, 0.3250, 0.0980];
c3=[0.9290, 0.6940, 0.1250];

figure
scatter(AvgMass1,AvgRProj1-0.005,200);
hold on;
scatter(AvgMass2,AvgRProj2+0.005,200);
scatter(AvgMass3,AvgRProj3,200);
Mx=0:0.2:20;
plot(Mx,3.3*exp(-(10/2/5.0)*(Mx/5).^1),'linewidth',2)
set(gca,'YScale','log')

text(AvgMass1-0.25,AvgRProj1-0.005,num2str(NumBin1),'Color',c1,'FontSize',14)
text(AvgMass2-0.25,AvgRProj2+0.005,num2str(NumBin2),'Color',c2,'FontSize',14)
text(AvgMass3-0.25,AvgRProj3+0.00,num2str(NumBin3),'Color',c3,'FontSize',14)

%axis([0.0 max(max(SimsAvgMass),max(DatAvgMass))*1.2 0 4])
xlabel('M bins [$10^5~M_\odot$]','interpreter','latex');
ylabel('$\left<r_\perp\right>$ [kpc]','interpreter','latex')
grid on;
set(gca,'FontSize',15);
%legend('Data','Simple dynamical friction model')
legend('Data Binning 1','Data Binning 2','Data Binning 3','Simple dynamical friction model')
%title(titleFig)
axis([0.5 16 0.1 4])
saveas(gcf,strcat('Plots/IntroCleanDataWithFit2.eps'),'epsc');

%% TRY SOME DISTANCE ERROR
DistError=1.2;
figure
scatter(AvgMass1*DistError^2,AvgRProj1*DistError-0.005,sizePoint1,'filled');
hold on;
scatter(AvgMass2*DistError^2,AvgRProj2*DistError+0.005,sizePoint2,'filled');
scatter(AvgMass3*DistError^2,AvgRProj3*DistError,sizePoint3,'filled');
Mx=0:0.2:20;
plot(Mx,3.6*exp(-(10/2/7)*(Mx/5).^1),'linewidth',2)
set(gca,'YScale','log')
%axis([0.0 max(max(SimsAvgMass),max(DatAvgMass))*1.2 0 4])
xlabel('M bins [$10^5~M_\odot$]','interpreter','latex');
ylabel('$\left<r_\perp\right>$ [kpc]','interpreter','latex')
grid on;
set(gca,'FontSize',15);
%legend('Data','Simple dynamical friction model')
legend('Data Binning 1','Data Binning 2','Data Binning 3','Simple dynamical friction model')
%title(titleFig)
axis([0.5 16 0.1 4])
saveas(gcf,strcat('Plots/IntroCleanDataWithFit3.eps'),'epsc');


%%
function [AvgMass,AvgRProj,sizePoint,NumPerBin]=getMomentsMR(MassBins)
BinDat = ReadDataIntoMassBins(MassBins);

binNum  = length(BinDat);
AvgMass = zeros(binNum,1);
AvgRProj = zeros(binNum,1);
sizePoint = zeros(binNum,1);
NumPerBin = zeros(binNum,1);
for ii = 1:binNum
    AvgMass(ii) = mean(BinDat{ii}(:,2));
    %ErrMass(ii) = std(BinDat{ii}(:,1));
    AvgRProj(ii) = mean(BinDat{ii}(:,1));
    %AvgR2Proj(ii) = mean((BinDat{ii}(:,1)).^2);
    sizePoint(ii) = 50*log(length(BinDat{ii}(:,1))/0.5)^1.5;
    NumPerBin(ii) = length(BinDat{ii}(:,1));
end

end