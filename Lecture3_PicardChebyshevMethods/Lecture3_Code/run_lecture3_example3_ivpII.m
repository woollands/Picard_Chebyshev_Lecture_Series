% Robyn Woollands 2017
% Texas A&M University - Department of Aerospace Engineering
% File name     : lecture3_example3_ivpII.m
% Description   : Second order Picard-Chebyshev IVP numerical integration for Lecture 3
%                 of the JUNKINS & WOOLLANDS, Picard-Chebyshev Lecture Series.
% Date Written  : April 1, 2017
% Date Modified : April 1, 2017
%================================================================

clear 
close all
clc

%% Integrate Second Order System: TPBVP (Example 3: Two-body Problem)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
mu      = 398600.4418;          % Gravitational Parameter (km^3 / s^2)
Re      = 6378.137;             % Earth Radius (km)
DU      = Re;                   % Canonical Distance Unit
TU      = sqrt(Re^3 / mu);      % Canonical Time Unit
VU      = DU/TU;                % Canonical Velocity Unit

% User can vary the following to generate different orbits...

N       = 40;                   % Chebyshev polynomial order
M       = N;                    % Number of sample points
tol     = 1e-13;                % Tolerance

a       = 8000;                 % Semimajor Axis (km)
Period  = 2*pi*sqrt(a^3/mu);    % Period (s)
e       = 0.125;                % Eccentricity
inc     = 10*pi/180;            % Inclination (rad)
w       = 0*pi/180;             % Argument of Perigee (rad)
Om      = 0*pi/180;             % Right Ascension of Ascending Node (rad)
p       = a*(1 - e^2);          % Semilatus Rectum (km)
MA      = 0*pi/180;             % Mean Anomaly (rad)
s       = 0*pi/180;             % "Special Parameter" (see Vallado)

t0      = 0;                    % Given initial time (s)
tf      = Period/4;             % Given final time (s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute r & v from elements
M0      = 0;
[r0,v0] = elm2rv(a,e,inc,Om,w,M0,s,mu); % "Given" initial conditions (km & km/s)

% Compute time vector
tau     = -cos((0:M)*pi/M);     % Cosine Sample Points
w1      = (tf + t0)/2;          % Scaling parameter
w2      = (tf - t0)/2;          % Scaling parameter
time    = (w2*tau + w1);        % Time (s)

%% Picard-Chebyshev
% Preparation
V0   = [v0(1) v0(2) v0(3); zeros(N-1,3)];                               % Initial condition vector (velocity)
X0   = [r0(1) r0(2) r0(3); zeros(N,3)];                                 % Initial condition vector (position)
X    = [r0(1).*ones(M+1,1)  r0(2).*ones(M+1,1) r0(3).*ones(M+1,1)];     % Cold start (initial position guess)
V    = [v0(1).*ones(M+1,1)  v0(2).*ones(M+1,1) v0(3).*ones(M+1,1)];     % Cold start (initial velocity guess)

% Scale factor (t -> tau)
w1   = (tf + t0)/2;
w2   = (tf - t0)/2;
% Time
time = w1 + w2*tau';

disp('%%% Picard-Chebyshev %%%')
tic
% Generate Picard-Chebysev constant matrices
[T2,P2,T1,P1,Ta,A] = clenshaw_curtis_ivpII(N,M);

itr  = 0;
err  = 10;
% Picard Integration (using Clenshaw & Curtis quadrature / path approximation)
while err > tol
        
    R3 = (X(:,1).^2 + X(:,2).^2 + X(:,3).^2).^(3/2);
    
    for i = 1:length(time)
        G(i,:) = -mu./R3(i).*X(i,:);    % Forcing function
    end
        
    % Velocity
    beta = w2.*P1*A*G + V0;             % Integration to obtain velocity coefficients
    Vnew = T1*beta;                     % Compute velocity from coefficients
    
    % Position
    alpha = w2.*P2*beta + X0;           % Integration to obtain position coefficients
    Xnew  = T2*alpha;                   % Compute position from coefficients
    
    % Compute error
    errX = max(max(abs(Xnew./DU - X./DU)));    % Position error
    errV = max(max(abs(Vnew./VU - V./VU)));    % Velocity error
    err  = max([errX errV]);                            % Combined error
    
    % Update
    V = Vnew;
    X = Xnew;
    
    % Iteration counter
    itr = itr + 1;
    if itr > 20
        disp(['Converged to: ',num2str(err),' instead of ',num2str(tol),'!!'])
        break
    end
    
end
toc

%% Analytical Solution (Truth)
for i = 1:length(time)
    [r, v] = FnG(t0,time(i),r0,v0,mu);
    FG(i,1:3) = r';
    FG(i,4:6) = v';
end

%% Plots
% Trajectory
figure(1)
plot3(FG(:,1)./DU,FG(:,2)./DU,FG(:,3)./DU,'r.-')
hold on
plot3(X(:,1)./DU,X(:,2)./DU,X(:,3)./DU,'g-*')
grid on
set(gca, 'FontName', 'Helvetica','FontSize',16)
xlabel('X (Re)')
ylabel('Y (Re)')
zlabel('Z (Re)')
legend('F&G Solution','Picard-Chebyshev')

figure(2)
subplot 311
plot(time./Period,FG(:,1)./DU,'k-','Linewidth',2)
hold on
plot(time./Period,X(:,1)./DU,'b*','markerSize',10)
set(gca, 'FontName', 'Helvetica','FontSize',16)
xlabel('Time (Periods)')
ylabel('X(t)')
legend('Truth','Picard-Chebyshev')
title('Position Components')

subplot 312
plot(time./Period,FG(:,2)./DU,'k-','Linewidth',2)
hold on
plot(time./Period,X(:,2)./DU,'b*','markerSize',10)
set(gca, 'FontName', 'Helvetica','FontSize',16)
xlabel('Time (Periods)')
ylabel('Y(t)')

subplot 313
plot(time./Period,FG(:,3)./DU,'k-','Linewidth',2)
hold on
plot(time./Period,X(:,3)./DU,'b*','markerSize',10)
set(gca, 'FontName', 'Helvetica','FontSize',16)
xlabel('Time (Periods)')
ylabel('Z(t)')

% Position Error
figure(3)
subplot 311
semilogy(time./Period,abs(FG(:,1)-X(:,1))./DU,'b*-','Linewidth',2)
set(gca, 'FontName', 'Helvetica','FontSize',16)
hold on
grid on
xlabel('Time (Periods)')
ylabel('|X - X_{truth}|/DU')
title('Position Error')

subplot 312
semilogy(time./Period,abs(FG(:,2)-X(:,2))./DU,'b*-','Linewidth',2)
set(gca, 'FontName', 'Helvetica','FontSize',16)
hold on
grid on
xlabel('Time (Periods)')
ylabel('|Y - Y_{truth}|/DU')

subplot 313
semilogy(time./Period,abs(FG(:,3)-X(:,3))./DU,'b*-','Linewidth',2)
set(gca, 'FontName', 'Helvetica','FontSize',16)
hold on
grid on
xlabel('Time (Periods)')
ylabel('|Z - Z_{truth}|/DU')

% Velocity Error
figure(4)
subplot 311
semilogy(time./Period,abs(FG(:,4)-V(:,1))./VU,'b*-','Linewidth',2)
set(gca, 'FontName', 'Helvetica','FontSize',16)
hold on
grid on
xlabel('Time (Periods)')
ylabel('|Vx - Vx_{truth}|/VU')
title('Velocity Error')

subplot 312
semilogy(time./Period,abs(FG(:,5)-V(:,2))./VU,'b*-','Linewidth',2)
set(gca, 'FontName', 'Helvetica','FontSize',16)
hold on
grid on
xlabel('Time (Periods)')
ylabel('|Vy - Vy_{truth}|/VU')

subplot 313
semilogy(time./Period,abs(FG(:,6)-V(:,3))./VU,'b*-','Linewidth',2)
set(gca, 'FontName', 'Helvetica','FontSize',16)
hold on
grid on
xlabel('Time (Periods)')
ylabel('|Vz - Vz_{truth}|/VU')
