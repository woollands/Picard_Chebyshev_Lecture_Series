% Robyn Woollands 2017
% Texas A&M University - Department of Aerospace Engineering
% File name     : lecture3_example1b_fvpI.m
% Description   : First order Picard-Chebyshev FVP numerical integration for Lecture 3 
%                 of the JUNKINS & WOOLLANDS, Picard-Chebyshev Lecture Series.
% Date Written  : March 26, 2017
% Date Modified : April 12, 2017
%================================================================

clear 
close all
clc

global eps
%% Integrate First Order System: IVP (Example 1: x_dot = -eps*x)

N    = 10;                  % Polynoimal order  
M    = 15;                  % Number of sample points
tau  = -cos([0:M].*pi/M);   % Cosine sample points
tol  = 1e-13;               % Tolerance

% Given
eps  = 0.05;                % Given system dynamics
xf   = 1;                   % Given initial conditions
t0   = 0;                   % Given initial time (s)
tf   = 25;                  % Given final time (s)

xf = 1;
%% Picard-Chebyshev
% Preparation
Xf           = [xf; zeros(N,1)];        % Initial condition vector
X            = xf.*ones(M+1,1);         % Initial solution guess (very poor guess still converge)

% Scale factor (t -> tau)
w1   = (tf + t0)/2;
w2   = -(tf - t0)/2;
% Time
time = w1 + w2*tau';

disp('%%% Picard-Chebyshev %%%')
tic
% Generate Picard-Chebysev constant matrices
[T1,P1,Ta,A] = clenshaw_curtis_fvpI(N,M); 

itr  = 0;
errX = 10;
% Picard Integration (using Clenshaw & Curtis quadrature / path approximation)
while errX > tol
    
    g    = -w2*eps.*X;      % Forcing function
    beta = P1*A*g + Xf;     % Integration to obtain coefficients
    Xnew = T1*beta;         % Compute solution from coefficients
    
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
x0 = X(end);
X_true = x0.*exp(-eps*time);

%% ode45
disp('%%% MATLAB: ode45 %%%')
tic
options = odeset('RelTol',tol,'AbsTol',tol);
[t,X_ode45] = ode45(@func_lecture3_example1,time,xf,options);
toc

%% Plots
% Trajectory
figure(1)
hold on
plot(time,X_true,'k-','Linewidth',2)
plot(time,X,'b*','MarkerSize',10)
plot(time,X_ode45,'r^','MarkerSize',10)
set(gca, 'FontName', 'Helvetica','FontSize',16)
xlabel('Time (s)')
ylabel('X(t)')
title(['X(t) = X0*exp(-',num2str(eps),'t)'])
legend('Truth','Picard-Chebyshev','ode45')

% Error
figure(2)
semilogy(time,abs(X_true - X),'b*-','Linewidth',2)
hold on
semilogy(time,abs(X_true - X_ode45),'ro-','Linewidth',2)
set(gca, 'FontName', 'Helvetica','FontSize',16)
grid on
xlabel('Time (s)')
ylabel('|X(t)_{true} - X(t)|')
title(['xdot = -',num2str(eps),'x'])
legend('Picard-Chebyshev','MATLAB: ode45')


