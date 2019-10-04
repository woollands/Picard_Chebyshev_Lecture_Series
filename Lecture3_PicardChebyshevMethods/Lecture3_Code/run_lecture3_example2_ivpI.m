% Robyn Woollands 2017
% Texas A&M University - Department of Aerospace Engineering
% File name     : lecture3_example2_ivpI.m
% Description   : First order Picard-Chebyshev IVP numerical integration for Lecture 3 
%                 of the JUNKINS & WOOLLANDS, Picard-Chebyshev Lecture Series.
% Date Written  : March 26, 2017
% Date Modified : April 12, 2017
%================================================================

clear 
close all
clc

global eps
%% Integrate First Order System: IVP (Example 1: x_dot = -eps*x)

N    = 100;                 % Polynoimal order  
M    = N;                   % Number of sample points
tau  = -cos([0:M].*pi/M);   % Cosine sample points
tol  = 1e-13;               % Tolerance

% Given
eps  = 0.05;                % Given system dynamics
x0   = 1;                   % Given initial conditions
t0   = 0;                   % Given initial time (s)
tf   = 25;                  % Given final time (s)

%% Picard-Chebyshev
% Preparation
X0           = [x0; zeros(N,1)];        % Initial condition vector
X            = x0.*ones(N+1,1);         % Initial solution guess (very poor guess still converge)

% Scale factor (t -> tau)
w1   = (tf + t0)/2;
w2   = (tf - t0)/2;
% Time
time = w1 + w2*tau';

disp('%%% Picard-Chebyshev %%%')
tic
% Generate Picard-Chebysev constant matrices
[T1,P1,Ta,A] = clenshaw_curtis_ivpI(N,M);  

itr  = 0;
errX = 10;
% Picard Integration (using Clenshaw & Curtis quadrature / path approximation)
while errX > tol
    
    g    = w2*cos(time + eps.*X);   % Forcing function
    beta = P1*A*g + X0;             % Integration to obtain coefficients
    Xnew = T1*beta;                 % Compute solution from coefficients
    
    % Compute error
    errX = max(abs(Xnew - X));
    
    % Update
    X = Xnew;
    
    % Iteration counter
    itr = itr + 1;
    if itr > 20
        disp(['Converged to: ',num2str(errX),' instead of ',num2str(tol),'!!'])
        break
    end
    
end
toc

%% Analytical Solution (Truth)
% Calculate parameters on page 38 (Bai's PhD dissertation) to get X_true.
gamma  = (eps / (1 + sqrt(1 - eps^2)));
alpha  = (2*eps) / (1 - eps + sqrt(1 - eps^2));
beta   = (1 / (1+alpha))*(tan((eps*x0) / 2));
phi    = (1/2)*(1 - gamma*eps)*time;
sigma  = alpha*(sin(phi) + beta*cos(phi));
X_true = -gamma*time + (2/eps)*atan2((beta + sigma.*cos(phi)), (1 + sigma.*sin(phi)));

%% ode45
disp('%%% MATLAB: ode45 %%%')
tic
options = odeset('RelTol',tol,'AbsTol',tol);
[t,X_ode45] = ode45(@func_lecture3_example2,time,x0,options);
toc

%% Plots
% Trajectory
figure(1)
hold on
plot(time,X_true,'k-','Linewidth',2)
plot(time,X,'b*','MarkerSize',10)
plot(time,X_ode45,'ro','MarkerSize',10)
set(gca, 'FontName', 'Helvetica','FontSize',16)
xlabel('Time (s)')
ylabel('X(t)')
legend('Truth','Picard-Chebyshev','MATLAB: ode45')

% Error
figure(2)
semilogy(time,abs(X_true - X),'b*-','Linewidth',2)
hold on
semilogy(time,abs(X_true - X_ode45),'ro-','Linewidth',2)
set(gca, 'FontName', 'Helvetica','FontSize',16)
grid on
xlabel('Time (s)')
ylabel('|X(t)_{true} - X(t)|')
title(['X(t) = cos(t +',num2str(eps),'X)'])
legend('Picard-Chebyshev','MATLAB: ode45')


