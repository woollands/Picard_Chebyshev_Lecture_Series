% Robyn Woollands 2017
% Texas A&M University - Department of Aerospace Engineering
% File name     : run_lecture1.m
% Description   : Runs functions for Lecture 1 of the JUNKINS & WOOLLANDS,
%                 Picard-Chebyshev Lecture Series
% Date Written  : March 16, 2017
% Date Modified : March 26, 2017
%================================================================

clear
close all
clc

%% Generate Chebyshev Polynomials (FIGURE 1)
N   = 10;       % Chebyshev polynomial order
M   = 100;      % Number of sample points
arg = 1;        % Use recursive or trigonometric Chebyshev formulation

[T,tau] = chebyshev(N,M,arg);

% Plot first 5 Chebyshev polynomials
figure(1)
hold on
plot(tau,T(:,1),'r-','Linewidth',3)
plot(tau,T(:,2),'b-','Linewidth',3)
plot(tau,T(:,3),'g-','Linewidth',3)
plot(tau,T(:,4),'k-','Linewidth',3)
plot(tau,T(:,5),'m-','Linewidth',3)
set(gca, 'FontName', 'Helvetica','FontSize',16)
title('Chebyshev Polynomials (First Kind)')
xlabel('\tau')
ylabel('T_{k}(\tau)')
legend('T_0','T_1','T_2','T_3','T_4')

%% Generate Cosine Samples (FIGURE 2)
clear

Mvec = [1 2 3 4 5 8 10 13 15];
for i = 1:length(Mvec)
    M = Mvec(i);
    
    % Cosine Sample Points
    tau = -cos([0:M].*pi/M);
    
    % Uniform Sample Points
    x = linspace(-1,1,M+1);
    
    % Plot sample points
    figure(2)
    hold on
    if i == 1
        plot(x,(M+0.2).*ones(M+1,1),'b.','MarkerSize',25)
        plot(tau,(M-0.2).*ones(M+1,1),'r.','MarkerSize',25)
    end
    plot(tau,(M+0.2).*ones(M+1,1),'k-','Linewidth',2)
    plot(x,(M+0.2).*ones(M+1,1),'b.','MarkerSize',25)
    plot(x,(M-0.2).*ones(M+1,1),'k-','Linewidth',2)
    plot(tau,(M-0.2).*ones(M+1,1),'r.','MarkerSize',25)
    
    if i == length(Mvec)
        set(gca, 'FontName', 'Helvetica','FontSize',16)
        title('Cosine vs Uniform Sampling')
        xlabel('\xi')
        ylabel('M (Integer Value)')
        legend('Uniform Sampling','Cosine Sampling')
        axis([-1.1 1.1 0 19])
    end
    
end

%% Demonstrate Runge Effect (FIGURE 3)
clear

% Truth
x       = linspace(-1,1,100);
y_true  = 1./(1 + 25*x.^2);

% Chebyshev Fit
N   = 10;       % Chebyshev polynomial order
M   = 10;       % Number of sample points
arg = 2;        % Use recursive or trigonometric Chebyshev formulation

[T,tau] = chebyshev(N,M,arg);

y_cheby     = 1./(1 + 25*tau'.^2);
p           = polyfit(tau',y_cheby,N);
y_cheby_fit = polyval(p,x);

% Power Series Fit
x_pow     = linspace(-1,1,N+1);
y_pow     = 1./(1 + 25*x_pow'.^2);
p         = polyfit(x_pow',y_pow,N);
y_pow_fit = polyval(p,x);

% Plot truth and approximations
figure(3)
hold on
% Functions
plot(x,y_true,'k','Linewidth',2)
plot(x,y_cheby_fit,'r','Linewidth',2)
plot(x,y_pow_fit,'b','Linewidth',2)
% Sample Points
plot(tau,y_cheby,'r.','MarkerSize',20)
plot(x_pow,y_pow,'b.','MarkerSize',20)
set(gca, 'FontName', 'Helvetica','FontSize',16)
title('Gaussian Function Approximation')
xlabel('\tau')
ylabel('y')
legend('Truth','Chebyshev','Power Series','Cosine Nodes','Uniform Nodes')

%% Chebyshev Least Squares Approximation to Ugly Function (FIGURES 4 & 5)
clear

Nvec    = [5:5:60];
for k = 1:length(Nvec);
    N       = Nvec(k);              % Chebyshev polynomial order
    M       = 110;                  % Number of sample points
    tau     = -cos((0:M)*pi/M);     % Cosine Sample Points
    
    % Compute time vector
    t0      = -1;
    tf      = 1;
    time    = tau;
    
    f_ugly  = @(x) x./2 + ((1./10 + x).*sin(5.*x - 1))./(1 + x.^2 .* sin(x-0.5).^2);
    
    G       = f_ugly(tau');
    
    % Approximate Acceleration
    [A,T]    = lsq_chebyshev_fit(N,M);  % Least Squares Operator (A)
    coeff    = A*G;                     % Least Squares Coefficients
    G_approx = T*coeff;                 % Approximated Acceleration
    
    figure(4)
    plot(tau,G,'b-','LineWidth',2)
    xlabel('x')
    ylabel('f(x)')
    title('f(x) = x/2 + ((1/10 + x)sin(5x - 1))/(1 + x^2sin(x-0.5)^2)')
    
    % Normalized Residuals
    for i = 1:M+1
        Err(i) = norm(G(i,:) - G_approx(i,:))/norm(G(1,:));
    end
    
    maxErr(k) = max(abs(Err));
    
end

figure(5)
semilogy(Nvec,maxErr,'r-','LineWidth',2)
hold on
semilogy(Nvec,maxErr,'b.','MarkerSize',20)
grid on
set(gca, 'FontName', 'Helvetica','FontSize',16)
xlabel('N')
ylabel('||G - G_{approx}|| / ||G_1||')
title('Fit Accuracy (Ugly Function)')

%% Chebyshev Least Squares Approximation (Two-Body Acceleration)
clear

mu      = 398600.4418;          % Gravitational Parameter
Re      = 6378.137;             % Earth Radius
DU      = Re;                   % Canonical Distance Unit
TU      = sqrt(Re^3 / mu);      % Canonical Time Unit

N       = 75;                   % Chebyshev polynomial order (N = 65 and 100 in Lecture 1)
M       = 110;                  % Number of sample points
tau     = -cos((0:M)*pi/M);     % Cosine Sample Points
Deg     = 70;                   % Spherical Harmonic Gravity Degree and Order (Deg = 40 and 70 in Lecture 1)

a       = 7000;                 % Semimajor Axis
Period  = 2*pi*sqrt(a^3/mu);    % Period
e       = 0.0;                  % Eccentricity
inc     = 45*pi/180;            % Inclination
w       = 0*pi/180;             % Argument of Perigee
Om      = 0*pi/180;             % Right Ascension of Ascending Node
p       = a*(1 - e^2);          % Semilatus Rectum
MA      = 0*pi/180;             % Mean Anomaly
s       = 0*pi/180;             % "Special Parameter" (see Vallado)

% Compute r & v from elements
M0      = 0;
[r0,v0] = elm2rv(a,e,inc,Om,w,M0,s,mu);

% Compute time vector
t0      = 0;
tf      = Period/3;
w1      = (tf + t0)/2;
w2      = (tf - t0)/2;
time    = (w2*tau + w1);

% Compute F&G Solution
for i = 1:length(time)
    [r, v] = FnG(t0,time(i),r0,v0,mu);
    X(i,:) = r';
    V(i,:) = v';
end

% Compute Acceleration
% NOTE:
% 1. gravitysphericalharmonic.m is a built in MATLAB function.
% 2. Older version of MATLAB do not contain this function.
% 3. Our research group at Texas A&M uses a more efficient C version with a MATLB wrapper.
% 4. The C version was written by A. BANI YOUNES, B. MACOMBER and A. PROBE.
[Gx, Gy, Gz] = gravitysphericalharmonic(X.*1e3,Deg);
G            = [Gx, Gy, Gz]./1e3;

% Approximate Acceleration
[A,T]    = lsq_chebyshev_fit(N,M);  % Least Squares Operator (A)
coeff    = A*G;                     % Least Squares Coefficients
G_approx = T*coeff;                 % Approximated Acceleration

% Normalized Residuals
for i = 1:M+1
    Err(i) = norm(G(i,:) - G_approx(i,:))/norm(G(1,:));
end

% Plot
figure(6)
semilogy(time./Period,Err,'r-','LineWidth',2)
hold on
grid on
set(gca, 'FontName', 'Helvetica','FontSize',16)
xlabel('time/Period')
ylabel('||G - G_{approx}|| / ||G_1||')
title('Acceleration Fit Accuracy')
% legend('M=110, N = 65, Deg = 40','M=110, N = 65, Deg = 70','M=110, N = 100, Deg = 70')

% Plot Earth
figure(7)
hold on
plot3(X(:,1)./Re,X(:,2)./Re,X(:,3)./Re,'yo-')
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')
title('Trajectory')

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
figure(7)
surface(x,y,z,props);
set(gca,'color','black')
