% Robyn Woollands 2017
% Texas A&M University - Department of Aerospace Engineering
% File name     : lecture3_example4d_ivpII_fvpII.m
% Description   : Second order Picard-Chebyshev BVP numerical integration for Lecture 3
%                 of the JUNKINS & WOOLLANDS, Picard-Chebyshev Lecture Series.
% Date Written  : March 27, 2017
% Date Modified : April 14, 2017
%================================================================

clear 
close all
clc

global mu Re omega Deg
%% Integrate Second Order System (Cascade): BVP (Example 4: Perturbed Two-body Problem)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
mu      = 398600.4418;          % Gravitational Parameter (km^3 / s^2)
Re      = 6378.137;             % Earth Radius (km)
omega   = 7.2921151e-5;         % Earth Rotation Rate
DU      = Re;                   % Canonical Distance Unit
TU      = sqrt(Re^3 / mu);      % Canonical Time Unit
VU      = DU/TU;                % Canonical Velocity Unit

% User can vary the following to generate different orbits...

N       = 100;                  % Chebyshev polynomial order
M       = 110;                  % Number of sample points
Deg     = 40;                   % Spherical Harmonic Gravity Degree and Order
tol     = 1e-13;                % Tolerance

a       = 20000;                % Semimajor Axis (km)
Period  = 2*pi*sqrt(a^3/mu);    % Period (s)
e       = 0.1;                  % Eccentricity
inc     = 20*pi/180;            % Inclination (rad)
w       = 0*pi/180;             % Argument of Perigee (rad)
Om      = 0*pi/180;             % Right Ascension of Ascending Node (rad)
p       = a*(1 - e^2);          % Semilatus Rectum (km)
MA      = 0*pi/180;             % Mean Anomaly (rad)
s       = 0*pi/180;             % "Special Parameter" (see Vallado)

t0      = 0;                    % Given initial time (s)
tf      = 0.09*Period;          % Given final time (s)

ivpfvp  = 0;                    % ivpfvp = 1 (v0 & rf), ivpfvp = 0 (vf & r0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute r & v from elements
M0      = 0;
[r0,v0] = elm2rv(a,e,inc,Om,w,M0,s,mu); % "Given" initial conditions (km & km/s)

% Compute time vector
tau     = -cos((0:M)*pi/M);     % Cosine Sample Points
w1      = (tf + t0)/2;          % Scaling parameter
w2_i    = (tf - t0)/2;          % Scaling parameter
w2_f    = -(tf - t0)/2;         % Scaling parameter
time    = (w2_i*tau + w1)';       % Time (s)

% Compute Analytical F&G Solution (initial guess or "warm start" for solving perturbed problem)
for i = 1:length(time)
    [r, v] = FnG(t0,time(i),r0,v0,mu);
    X(i,:) = r';
    V(i,:) = v';
end
FG = [X V];

%% Picard-Chebyshev
% Preparation
if ivpfvp == 1
    V0   = [FG(1,4:6); zeros(N-1,3)];     % Initial condition vector (velocity)
    Xf   = [FG(end,1:3); zeros(N,3)];     % Final condition vector (position)
elseif ivpfvp == 0
    X0   = [FG(1,1:3); zeros(N,3)];       % Initial condition vector (position)
    Vf   = [FG(end,4:6); zeros(N-1,3)];   % Final condition vector (velocity)
end

disp('%%% Picard-Chebyshev %%%')
tic
% Generate Picard-Chebysev constant matrices
[T2_i,P2_i,T1_i,P1_i,Ta_i,A_i] = clenshaw_curtis_ivpII(N,M);
[T2_f,P2_f,T1_f,P1_f,Ta_f,A_f] = clenshaw_curtis_fvpII(N,M);

itr  = 0;
err = 10;vec = [];
% Picard Integration (using Clenshaw & Curtis quadrature / path approximation)
while err > tol
    
    % Convert to ECEF
    [xB,~] = eci2ecef(time,X,0,omega);
        
    for i = 1:length(time)
        Gout(i,:) = EGMGravMex([xB(i,1), xB(i,2), xB(i,3)].*1e3,Deg)./1e3;     % Forcing function
    end
    
    % Convert to Inertial
    G = ecef2eci(time,Gout,omega);

    % Given initial velocity and final position (IVP then FVP)
    if ivpfvp == 1;
        beta  = w2_i.*P1_i*A_i*G + V0;    % Velocity coefficients
        Vnew  = T1_i*beta;                % Velocity
        alpha = w2_f.*P2_f*beta + Xf;     % Position coefficients
        Xnew  = T2_i*alpha;               % Position
    end
 
    % Given final velocity and initial position (FVP then IVP)
    if ivpfvp == 0;
        beta  = w2_f.*P1_f*A_i*G + Vf;    % Velocity coefficients
        Vnew  = T1_i*beta;                        % Velocity
        alpha = w2_i.*P2_i*beta + X0;             % Position coefficients
        Xnew  = T2_i*alpha;                       % Position
    end
    
    % Compute error
    errX = max(max(abs(Xnew./DU - X./DU)));    % Position error
    errV = max(max(abs(Vnew./VU - V./VU)));    % Velocity error
    err  = max([errX errV]);                   % Combined error
    
    % Update
    V = Vnew;
    X = Xnew;
    
    % Iteration counter
    itr = itr + 1;
    if itr > 20
        disp(['Converged to: ',num2str(err),' instead of ',num2str(tol),'!!'])
        break
    end
    
    if ivpfvp == 1
        vec  = [vec;abs([Xnew(1,:) Vnew(end,:)])];
    elseif ivpfvp == 0
        vec  = [vec;abs([Xnew(end,:) Vnew(1,:)])];
    end
    
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
% axis equal
% legend('F&G Solution','Picard-Chebyshev')

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

% Trajectory Components and Boundary Conditions
if ivpfvp == 1
    figure(3)
    subplot 311
    hold on
    plot(time./Period,X(:,1)./DU,'b-','Linewidth',2)
    plot(time(end)./Period,Xf(1,1)./DU,'r.','MarkerSize',20)
    plot(time(1)./Period,X(1,1)./DU,'g.','MarkerSize',20)
    set(gca, 'FontName', 'Helvetica','FontSize',16)
    xlabel('Time (Periods)')
    ylabel('X (DU)')
%     legend('Trajectory','Specified Final BC','Converged Initial BC')
    subplot 312
    hold on
    plot(time./Period,X(:,2)./DU,'b-','Linewidth',2)
    plot(time(end)./Period,Xf(1,2)./DU,'r.','MarkerSize',20)
    plot(time(1)./Period,X(1,2)./DU,'g.','MarkerSize',20)
    set(gca, 'FontName', 'Helvetica','FontSize',16)
    xlabel('Time (Periods)')
    ylabel('Y (DU)')
%     legend('Trajectory','Specified Final BC','Converged Initial BC')
    subplot 313
    hold on
    plot(time./Period,X(:,3)./DU,'b-','Linewidth',2)
    plot(time(end)./Period,Xf(1,3)./DU,'r.','MarkerSize',20)
    plot(time(1)./Period,X(1,3)./DU,'g.','MarkerSize',20)
    set(gca, 'FontName', 'Helvetica','FontSize',16)
    xlabel('Time (Periods)')
    ylabel('Z (DU)')
%     legend('Trajectory','Specified Final BC','Converged Initial BC')
    
    figure(4)
    subplot 311
    hold on
    plot(time./Period,V(:,1)./VU,'b-','Linewidth',2)
    plot(time(1)./Period,V0(1,1)./VU,'r.','MarkerSize',20)
    plot(time(end)./Period,V(end,1)./VU,'g.','MarkerSize',20)
    set(gca, 'FontName', 'Helvetica','FontSize',16)
    xlabel('Time (Periods)')
    ylabel('Vx (DU/TU)')
%     legend('Trajectory','Specified Initial BC','Converged Final BC')
    subplot 312
    hold on
    plot(time./Period,V(:,2)./VU,'b-','Linewidth',2)
    plot(time(1)./Period,V0(1,2)./VU,'r.','MarkerSize',20)
    plot(time(end)./Period,V(end,2)./VU,'g.','MarkerSize',20)
    set(gca, 'FontName', 'Helvetica','FontSize',16)
    xlabel('Time (Periods)')
    ylabel('Vy (DU/TU)')
%     legend('Trajectory','Specified Initial BC','Converged Final BC')
    subplot 313
    hold on
    plot(time./Period,V(:,3)./VU,'b-','Linewidth',2)
    plot(time(1)./Period,V0(1,3)./VU,'r.','MarkerSize',20)
    plot(time(end)./Period,V(end,3)./VU,'g.','MarkerSize',20)
    set(gca, 'FontName', 'Helvetica','FontSize',16)
    xlabel('Time (Periods)')
    ylabel('Vz (DU/TU)')
%     legend('Trajectory','Specified Initial BC','Converged Final BC')
    
elseif ivpfvp == 0
    figure(3)
    subplot 311
    hold on
    plot(time./Period,X(:,1)./DU,'b-','Linewidth',2)
    plot(time(1)./Period,X0(1,1)./DU,'r.','MarkerSize',20)
    plot(time(end)./Period,X(end,1)./DU,'g.','MarkerSize',20)
    set(gca, 'FontName', 'Helvetica','FontSize',16)
    xlabel('Time (Periods)')
    ylabel('X (DU)')
%     legend('Trajectory','Specified Initial BC','Converged Final BC')
    subplot 312
    hold on
    plot(time./Period,X(:,2)./DU,'b-','Linewidth',2)
    plot(time(1)./Period,X0(1,2)./DU,'r.','MarkerSize',20)
    plot(time(end)./Period,X(end,2)./DU,'g.','MarkerSize',20)
    set(gca, 'FontName', 'Helvetica','FontSize',16)
    xlabel('Time (Periods)')
    ylabel('Y (DU)')
%     legend('Trajectory','Specified Initial BC','Converged Final BC')
    subplot 313
    hold on
    plot(time./Period,X(:,3)./DU,'b-','Linewidth',2)
    plot(time(1)./Period,X0(1,3)./DU,'r.','MarkerSize',20)
    plot(time(end)./Period,X(end,3)./DU,'g.','MarkerSize',20)
    set(gca, 'FontName', 'Helvetica','FontSize',16)
    xlabel('Time (Periods)')
    ylabel('Z (DU)')
%     legend('Trajectory','Specified Initial BC','Converged Final BC')
    
    figure(4)
    subplot 311
    hold on    
    plot(time./Period,V(:,1)./VU,'b-','Linewidth',2)
    plot(time(end)./Period,Vf(1,1)./VU,'r.','MarkerSize',20)
    plot(time(1)./Period,V(1,1)./VU,'g.','MarkerSize',20)
    set(gca, 'FontName', 'Helvetica','FontSize',16)
    xlabel('Time (Periods)')
    ylabel('Vx (DU/TU)')
%     legend('Trajectory','Specified Final BC','Converged Initial BC')
    subplot 312
    hold on
    plot(time./Period,V(:,2)./VU,'b-','Linewidth',2)
    plot(time(end)./Period,Vf(1,2)./VU,'r.','MarkerSize',20)
    plot(time(1)./Period,V(1,2)./VU,'g.','MarkerSize',20)
    set(gca, 'FontName', 'Helvetica','FontSize',16)
    xlabel('Time (Periods)')
    ylabel('Vy (DU/TU)')
%     legend('Trajectory','Specified Final BC','Converged Initial BC')
    subplot 313
    hold on
    plot(time./Period,V(:,3)./VU,'b-','Linewidth',2)
    plot(time(end)./Period,Vf(1,3)./VU,'r.','MarkerSize',20)
    plot(time(1)./Period,V(1,3)./VU,'g.','MarkerSize',20)
    set(gca, 'FontName', 'Helvetica','FontSize',16)
    xlabel('Time (Periods)')
    ylabel('Vz (DU/TU)')
%     legend('Trajectory','Specified Final BC','Converged Initial BC')
end

if ivpfvp == 1
   figure(5)
   subplot 311
   hold on
   plot([1:itr-1],(vec(:,1)-Xnew(1,1)).*1e3,'k-','Linewidth',2)
   plot(1,(vec(1,1)-Xnew(1,1)).*1e3,'m.','MarkerSize',20)
   plot(itr-1,0,'g.','MarkerSize',20)
   set(gca, 'FontName', 'Helvetica','FontSize',16)
   ylabel('X^i_0 - X^c_0 (m)')
   xlabel('Iterations')
%    legend('Convergence Path','Two-body Initial Guess X_0','Converged X_0')
   subplot 312
   hold on
   plot([1:itr-1],(vec(:,2)-Xnew(1,2)).*1e3,'k-','Linewidth',2)
   plot(1,(vec(1,2)-Xnew(1,2)).*1e3,'m.','MarkerSize',20)
   plot(itr-1,0,'g.','MarkerSize',20)
   set(gca, 'FontName', 'Helvetica','FontSize',16)
   ylabel('Y^i_0 - Y^c_0 (m)')
   xlabel('Iterations')
%    legend('Convergence Path','Two-body Initial Guess Y_0','Converged Y_0')
   subplot 313
   hold on
   plot([1:itr-1],(vec(:,3)-Xnew(1,3)).*1e3,'k-','Linewidth',2)
   plot(1,(vec(1,3)-Xnew(1,3)).*1e3,'m.','MarkerSize',20)
   plot(itr-1,0,'g.','MarkerSize',20)
   set(gca, 'FontName', 'Helvetica','FontSize',16)
   ylabel('Z^i_0 - Z^c_0 (m)')
   xlabel('Iterations')
%    legend('Convergence Path','Two-body Initial Guess Z_0','Converged Z_0')
   
   figure(6)
   subplot 311
   hold on
   plot([1:itr-1],(vec(:,4)+Vnew(end,1)).*1e3,'k-','Linewidth',2)
   plot(1,(vec(1,4)+Vnew(end,1)).*1e3,'m.','MarkerSize',20)
   plot(itr-1,0,'g.','MarkerSize',20)
   set(gca, 'FontName', 'Helvetica','FontSize',16)
   ylabel('Vx^i_f - Vx^c_f (m/s)')
   xlabel('Iterations')
%    legend('Convergence Path','Two-body Initial Guess Vx_0','Converged Vx_0')
   subplot 312
   hold on
   plot([1:itr-1],(vec(:,5)-Vnew(end,2)).*1e3,'k-','Linewidth',2)
   plot(1,(vec(1,5)-Vnew(end,2)).*1e3,'m.','MarkerSize',20)
   plot(itr-1,0,'g.','MarkerSize',20)
   set(gca, 'FontName', 'Helvetica','FontSize',16)
   ylabel('Vy^i_f - Vy^c_f (m/s)')
   xlabel('Iterations')
%    legend('Convergence Path','Two-body Initial Guess Vy_0','Converged Vy_0')
   subplot 313
   hold on
   plot([1:itr-1],(vec(:,6)-Vnew(end,3)).*1e3,'k-','Linewidth',2)
   plot(1,(vec(1,6)-Vnew(end,3)).*1e3,'m.','MarkerSize',20)
   plot(itr-1,0,'g.','MarkerSize',20)
   set(gca, 'FontName', 'Helvetica','FontSize',16)
   ylabel('Vz^i_f - Vz^c_f (m/s)')
   xlabel('Iterations')
%    legend('Convergence Path','Two-body Initial Guess Vz_0','Converged Vz_0')
   
elseif ivpfvp == 0
   figure(7)
   subplot 311
   hold on
   plot([1:itr],(vec(:,1)-Xnew(end,1)).*1e3,'k-','Linewidth',2)
   plot(1,(vec(1,1)-Xnew(end,1)).*1e3,'m.','MarkerSize',20)
   plot(itr,0,'g.','MarkerSize',20)
   set(gca, 'FontName', 'Helvetica','FontSize',16)
   ylabel('X^i_f - X^c_f (m)')
   xlabel('Iterations')
%    legend('Convergence Path','Two-body Initial Guess X_f','Converged X_f')
   subplot 312
   hold on
   plot([1:itr],(vec(:,2)-Xnew(end,2)).*1e3,'k-','Linewidth',2)
   plot(1,(vec(1,2)-Xnew(end,2)).*1e3,'m.','MarkerSize',20)
   plot(itr,0,'g.','MarkerSize',20)
   set(gca, 'FontName', 'Helvetica','FontSize',16)
   ylabel('Y^i_f - Y^c_f (m)')
   xlabel('Iterations')
%    legend('Convergence Path','Two-body Initial Guess Y_f','Converged Y_f')
   subplot 313
   hold on
   plot([1:itr],(vec(:,3)-Xnew(end,3)).*1e3,'k-','Linewidth',2)
   plot(1,(vec(1,3)-Xnew(end,3)).*1e3,'m.','MarkerSize',20)
   plot(itr,0,'g.','MarkerSize',20)
   set(gca, 'FontName', 'Helvetica','FontSize',16)
   ylabel('Z^i_0 - Z^c_0 (m)')
   xlabel('Iterations')
%    legend('Convergence Path','Two-body Initial Guess Z_f','Converged Z_f')
   
   figure(8)
   subplot 311
   hold on
   plot([1:itr],(vec(:,4)-Vnew(1,1)).*1e3,'k-','Linewidth',2)
   plot(1,(vec(1,4)-Vnew(1,1)).*1e3,'m.','MarkerSize',20)
   plot(itr,0,'g.','MarkerSize',20)
   set(gca, 'FontName', 'Helvetica','FontSize',16)
   ylabel('Vx^i_0 - Vx^c_0 (m/s)')
   xlabel('Iterations')
%    legend('Convergence Path','Two-body Initial Guess Vx_f','Converged Vx_f')
   subplot 312
   hold on
   plot([1:itr],(vec(:,5)-Vnew(1,2)).*1e3,'k-','Linewidth',2)
   plot(1,(vec(1,5)-Vnew(1,2)).*1e3,'m.','MarkerSize',20)
   plot(itr,0,'g.','MarkerSize',20)
   set(gca, 'FontName', 'Helvetica','FontSize',16)
   ylabel('Vy^i_f - Vy^c_f (m/s)')
   xlabel('Iterations')
%    legend('Convergence Path','Two-body Initial Guess Vy_f','Converged Vy_f')
   subplot 313
   hold on
   plot([1:itr],(vec(:,6)-Vnew(1,3)).*1e3,'k-','Linewidth',2)
   plot(1,(vec(1,6)-Vnew(1,3)).*1e3,'m.','MarkerSize',20)
   plot(itr,0,'g.','MarkerSize',20)
   set(gca, 'FontName', 'Helvetica','FontSize',16)
   ylabel('Vz^i_f - Vz^c_f (m/s)')
   xlabel('Iterations')
%    legend('Convergence Path','Two-body Initial Guess Vz_f','Converged Vz_f')
end
