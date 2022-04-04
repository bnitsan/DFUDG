function [dx] = nbodyCSoft(t,x)
global G Gn m rS densityS massS sigmaS TauFudgeFac epsilon LogOption

N = length(m); % The number of bodies

% computing distances r between bodies
X=x(1:N);
Y=x(N+1:2*N);
Z=x(2*N+1:3*N);
VX=x(3*N+1:4*N);
VY=x(4*N+1:5*N);
VZ=x(5*N+1:6*N);

DX=distance(X); % external function below
DY=distance(Y);
DZ=distance(Z);

%r  = (DX.^2+DY.^2+DZ.^2).^0.5; % matrix of distances between all particles
rSoft  = (DX.^2+DY.^2+DZ.^2+epsilon^2).^0.5; % matrix of distances between all particles WITH SOFTENING
rc = (X.^2+Y.^2+Z.^2).^0.5;    % vector of distances from center
V  = (VX.^2+VY.^2+VZ.^2).^0.5;    % vector of velocities
sigma=interp1(rS,sigmaS,rc);
Xv  = V./sigma/sqrt(2);

rho  = interp1(rS,densityS,rc);
MEnc = interp1(rS,massS,rc);
TauFudgeFacLocal = interp1(rS,TauFudgeFac,rc);
if LogOption == 1 
    Taus        = 0.789.*V.^3./(m.*rho.*...
        (erf(Xv)-2*Xv.*exp(-Xv.^2)/sqrt(pi)).*...
        log(1+(2*rc.*V.^2./m/0.449).^2))+TauFudgeFacLocal; % ISO-like log
else
    Taus        = 0.789.*V.^3./(m.*rho.*...
        (erf(Xv)-2*Xv.*exp(-Xv.^2)/sqrt(pi)).*...
        log(1+(0.5*sigma.^2./m/0.449).^2))+TauFudgeFacLocal; % NFW-like log
end

% accelerations
M=ones(N,1)*m'; % square matrix of masses, each row is m(i) repeated along columns (for next line)
%M = M';
%M;
%m*ones(1,N);
a=-Gn*M./rSoft.^2; % accelerations of each mass
   
% get rid of diagonal entries (Inf - forces on bodies on itself)
a(1:N+1:N*N) = 0;
    
ax =a.* (DX./ rSoft);
ay =a.* (DY./ rSoft);
az =a.* (DZ./ rSoft);
   
ax(1:N+1:N*N) = 0;	
ay(1:N+1:N*N) = 0;
az(1:N+1:N*N) = 0;
   
% Mean field 
   
a_xC=-G*MEnc.*X./rc.^3;
a_yC=-G*MEnc.*Y./rc.^3;
a_zC=-G*MEnc.*Z./rc.^3;
 
% Dynamical friction
   
a_DF_x = -VX./Taus;
a_DF_y = -VY./Taus;
a_DF_z = -VZ./Taus;

% net accel on mass m(i) -- adding all a's along column number (direction 2)
a_sx=sum(ax,2)'+a_xC'+a_DF_x';
a_sy=sum(ay,2)'+a_yC'+a_DF_y';
a_sz=sum(az,2)'+a_zC'+a_DF_z';
   
dx = [x(3*N+1:6*N)' a_sx a_sy a_sz]';

%size(Taus)
%size(dx)
end

function [dist]=distance(X)
% this function calculates distance between two coordinates in a given vector X
dist = X-X.'; 
end
