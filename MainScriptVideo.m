modelName = 'NFWMLS';modelNum  = '1'; saveName=strcat(modelName,modelNum,'s1A');
load(strcat('Sims',modelName,modelNum,'.mat'))

%%
mInitial=St1A{1}.mtab(1,:)';
log(mInitial)/max(log(mInitial))
%%
plotVideo(St1A{1},true,'vid',5,2,1.0)
%%
plotVideo(St1A{2},false,'vid',5,2,1.0)

%%
function plotVideo(sim,recordFlag,VidName,maxR,viewType,ColorFactor)
TTT = sim.TTT;
XXX = sim.XXX;
mtab = sim.mtab;
N=(size(XXX,2)/6);
i=1:N;
XX = []; YY = []; ZZ = [];
%VX = []; VY = []; VZ = [];
XX(:,i)=XXX(:,i);
YY(:,i)=XXX(:,i+N);
ZZ(:,i)=XXX(:,i+2*N);

if recordFlag
    vid = VideoWriter(strcat('Plots/',VidName,'.mp4'),'MPEG-4');
    open(vid);
end
myplot=figure('Position',[100 100 850 850]);hold on; % figure window
set(gca,'FontSize',15);
axis([-maxR maxR -maxR maxR -maxR maxR]*1);
hold on;
xlabel('X [kpc]');ylabel('Y [kpc]');zlabel('Z [kpc]');
daspect([1 1 1]);
switch viewType
    case 1
        view(0,0)
        time  = text(-maxR*2/3,0,-maxR*2/3,['t [Gyr] = ' num2str(0)],'FontSize',14);
        GCNum = text(-maxR*2/3,0,1-maxR*2.5/3,['#GC = ' num2str(N)],'FontSize',14);
    case 2
        view(0, 90)
        time  = text(-maxR*2/3,-maxR*2/3,0,['t [Gyr] = ' num2str(0)],'FontSize',22);
        GCNum = text(-maxR*2/3,-maxR*2.5/3,0,['#GC = ' num2str(N)],'FontSize',22);
end
%view(2); % change to view(3) for 3D projection
grid on
%whitebg('black')

mInitial=mtab(1,:)';
interpMass= 0.9*mInitial/max(mInitial);

Ncolors=20;
colorPos=parula(Ncolors);
colorMass=ColorFactor*interp1(0:(1/(Ncolors-1)):1,colorPos,interpMass);
colorbar
colormap parula
% plots for each orbit of each mass
X = XX(1,:); Y = YY(1,:); Z = ZZ(1,:);
for ii=1:N
   body_traj(ii)=plot3(X(ii), Y(ii),Z(ii),'-','Color',colorMass(ii,:)); hold on       % orbit line
   body(ii)=plot3(X(ii), Y(ii), Z(ii),'o', 'MarkerSize',15,'Color',colorMass(ii,:),'MarkerFaceColor',colorMass(ii,:));    % mass marker
end
drawnow

if recordFlag
    A=getframe(myplot);
    writeVideo(vid,A);
end

deltaT=0.01;
NT=round(TTT(end)/deltaT);
TTTind=zeros(length(NT),1);
TTTind(1)=1;
for k=2:NT
    TTTind(k)=find(TTT>=deltaT*(k-1),1,'first');
end
TTTind(length(TTTind)+1) = length(TTT);

for k=2:length(TTTind)
    time.String=['t [Gyr] = ' num2str(round(TTT(TTTind(k)),3))];
    GCNum.String=['#GC = ' num2str(sum(mtab(TTTind(k),:)>1e-8))];

    mCurrent=mtab(TTTind(k),:)';
    interpMass= mCurrent/max(mCurrent)/1.1;
    colorMass=ColorFactor*interp1(0:(1/(Ncolors-1)):1,colorPos,interpMass);
    for ii=1:N
        if mCurrent(ii)>1e-8
            %disp(colorMass)
            set(body(ii), 'XData',XX(TTTind(k),ii),'YData',YY(TTTind(k),ii),'ZData',ZZ(TTTind(k),ii),'MarkerSize',15,'Color',colorMass(ii,:)); % update mass' markers
        else
            set(body(ii), 'XData',100,'YData',100,'ZData',ZZ(TTTind(k),ii),'Color',colorMass(ii,:)); % update mass' markers
        end    
    end
    pause(0.002)
    
    if recordFlag
        A=getframe(myplot);
        writeVideo(vid,A);
    end

    %drawnow
end

if recordFlag
    %A=getframe(myplot);
    %writeVideo(vid,A);
    close(vid); % close video file
end
hold off
end
