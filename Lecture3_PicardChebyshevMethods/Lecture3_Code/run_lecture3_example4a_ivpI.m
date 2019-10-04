% Robyn Woollands 2017
% Texas A&M University - Department of Aerospace Engineering
% File name     : lecture3_example4a_ivpI.m
% Description   : First order Picard-Chebyshev IVP numerical integration for Lecture 3
%                 of the JUNKINS & WOOLLANDS, Picard-Chebyshev Lecture Series.
% Date Written  : March 27, 2017
% Date Modified : April 12, 2017
%================================================================

clear 
close all
clc

global mu Re omega Deg
%% Integrate Second Order System in First Order Form: IVP (Example 4: Perturbed Two-body Problem)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
mu      = 398600.4418;          % Gravitational Parameter (km^3 / s^2)
Re      = 6378.137;             % Earth Radius (km)
omega   = 7.2921151e-5;         % Earth Rotation Rate
DU      = Re;                   % Canonical Distance Unit
TU      = sqrt(Re^3 / mu);      % Canonical Time Unit
VU      = DU/TU;                % Canonical Velocity Unit

% User can vary the following to generate different orbits...

N       = 80;                  % Chebyshev polynomial order
M       = N;                  % Number of sample points
Deg     = 40;                   % Spherical Harmonic Gravity Degree and Order
tol     = 1e-13;                % Tolerance

a       = 8000;                 % Semimajor Axis (km)
Period  = 2*pi*sqrt(a^3/mu);    % Period (s)
e       = 0.125;                % Eccentricity
inc     = 20*pi/180;            % Inclination (rad)
w       = 0*pi/180;             % Argument of Perigee (rad)
Om      = 0*pi/180;             % Right Ascension of Ascending Node (rad)
p       = a*(1 - e^2);          % Semilatus Rectum (km)
MA      = 0*pi/180;             % Mean Anomaly (rad)
s       = 0*pi/180;             % "Special Parameter" (see Vallado)

t0      = 0;                    % Given initial time (s)
tf      = Period/3;             % Given final time (s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute r & v from elements
M0      = 0;
[r0,v0] = elm2rv(a,e,inc,Om,w,M0,s,mu); % "Given" initial conditions (km & km/s)

% Compute time vector
tau     = -cos((0:M)*pi/M);     % Cosine Sample Points
w1      = (tf + t0)/2;          % Scaling parameter
w2      = (tf - t0)/2;          % Scaling parameter
time    = (w2*tau + w1)';       % Time (s)

% Compute Analytical F&G Solution (initial guess or "warm start" for solving perturbed problem)
for i = 1:length(time)
    [r, v] = FnG(t0,time(i),r0,v0,mu);
    X(i,:) = r';
    V(i,:) = v';
end
FG = [X V];

%% Picard-Chebyshev
% Preparation
IC   = [[r0' v0']; zeros(N,6)];   % Initial condition vector (velocity)

disp('%%% Picard-Chebyshev %%%')
tic
% Generate Picard-Chebysev constant matrices
[T1,P1,Ta,A] = clenshaw_curtis_ivpI(N,M);

itr    = 0;
err    = 10;
errVec = [];
% Picard Integration (using Clenshaw & Curtis quadrature / path approximation)
while err > tol
    
    % Convert to ECEF
    [xB,~] = eci2ecef(time,X,0,omega);
        
    % Compute Acceleration
    for i = 1:length(time)
        Gout(i,:) = EGMGravMex([xB(i,1), xB(i,2), xB(i,3)].*1e3,Deg)./1e3;
    end
    
    % Convert to Inertial
    G = ecef2eci(time,Gout,omega);
        
    % Delta States
    F = [V G];
    
    % States
    beta   = w2.*P1*A*F + IC;   % Coefficients of states 
    states = T1*beta;           % States
    
    % Compute error
    errX = max(max(abs(states(:,1:3)./DU - X./DU)));    % Position error
    errV = max(max(abs(states(:,4:6)./VU - V./VU)));    % Velocity error
    err  = max([errX errV]);                            % Combined error                      
    
    % Update
    X = states(:,1:3);
    V = states(:,4:6);
    
    % Iteration counter
    itr = itr + 1;
    if itr > 30
        disp(['Converged to: ',num2str(err),' instead of ',num2str(tol),'!!'])
        break
    end
    
    % Store Error
    errVec = [errVec; [errX errV err]];
    
end
toc

% Convert to ECEF
[xB,vB] = eci2ecef(time,X,V,omega);

% Compute Hamiltonian
H       = jacobi_integral([xB vB].*1e3);

%% Plots
figure(1)
plot3(FG(:,1)./DU,FG(:,2)./DU,FG(:,3)./DU,'r.-')
hold on
plot3(X(:,1)./DU,X(:,2)./DU,X(:,3)./DU,'y-o')
grid on
set(gca, 'FontName', 'Helvetica','FontSize',16)
xlabel('X (Re)')
ylabel('Y (Re)')
zlabel('Z (Re)')
legend('F&G Solution','Picard-Chebyshev')

% NOTE: The following code (17 lines) was obtained from MathWorks
% online file exchange (Ryan Gray).

load('topo.mat','topo','topomap1');
whos topo topomap1;
colormap(topomap1);
% Create the surface.
[x,y,z] = sphere(50);
props.AmbientStrength           = 0.1;
props.DiffuseStrength           = 1;
props.SpecularColorReflectance  = .5;
props.SpecularExponent          = 20;
props.SpecularStrength          = 1;
props.FaceColor                 = 'texture';
props.EdgeColor                 = 'none';
props.FaceLighting              = 'phong';
props.Cdata                     = topo;
figure(1)
surface(x,y,z,props);
set(gca,'color','black')

% Hamiltonian 
figure(2)
semilogy(time./Period,H,'r.-','linewidth',2);
grid on
set(gca, 'FontName', 'Helvetica','FontSize',16)
title('Hamiltonian')
xlabel('Time (Periods)')
ylabel('| E - E_o | / | E_o |')

% Iterations
figure(3)
semilogy([1:itr],errVec(:,3),'k.-','Linewidth',2,'MarkerSize',20)
hold on
grid on
semilogy([1:itr],errVec(:,1),'b-','Linewidth',2)
semilogy([1:itr],errVec(:,2),'r-','Linewidth',2)
set(gca, 'FontName', 'Helvetica','FontSize',16)
xlabel('Iterations')
ylabel('Relative Convergence Error')
legend('Combined Relative Error','Position Relative Error','Velocity Relative Error')

