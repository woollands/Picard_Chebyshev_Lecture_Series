% Robyn Woollands 2017
% Texas A&M University - Department of Aerospace Engineering
% File name     : run_lecture2.m
% Description   : Runs functions for Lecture 2 of the JUNKINS & WOOLLANDS,
%                 Picard-Chebyshev Lecture Series
% Date Written  : March 18, 2017
% Date Modified : March 18, 2017
%================================================================

clear 
close all
clc

% Limits of integration;
a          = -1;  % Lower 
b          = 1;   % Upper
% Sample points (for plotting)
x          = linspace(a,b,100); % Uniform sampling
% Desired tolerance
tol        = 1e-10;
% Functions
f_ugly     = @(x) x./2 + ((1./10 + x).*sin(5.*x - 1))./(1 + x.^2 .* sin(x-0.5).^2);
f_exp      = @(x) (4.*x + 4).*exp(4.*x + 4);
% Truth
truth_exp  = 5216.926477;
truth_ugly = 0.022619973769678;

%% Quadrature (MATLAB quad)
[q_ugly,fevals_ugly] = quad(f_ugly,a,b,tol);
[q_exp,fevals_exp]   = quad(f_exp,a,b,tol);

%% Gaussian Quadrature (Hard coded)
method = [2 3 4];
for i = 1:3
    I_exp(i)  = gaussian_quad(method(i),f_exp);
    I_ugly(i) = gaussian_quad(method(i),f_ugly);
end

%% Clenshaw & Curtis Quadrature
N            = 20;                          % Polynoimal order  
M            = 20;                           % Number of sample points
tau          = -cos([0:M].*pi/M);           % Cosine sample points
[T1,P1,Ta,A] = clenshaw_curtis_ivpI(N,M);   % Constant matrices
% Ugly function
beta_ugly    = P1*A*f_ugly(tau');           % Coefficients
X_ugly       = T1*beta_ugly;                % Path approximation
I_CN_ugly    = X_ugly(end);                 % Quadraeture
% Exponential function
beta_exp     = P1*A*f_exp(tau');            % Coefficients
X_exp        = T1*beta_exp;                 % Path approximation
I_CN_exp     = X_exp(end);                  % Quadraeture

%% Compute Errors
% Errors Ugly
ErrU(1) = abs(truth_ugly-q_ugly);
ErrU(2) = abs(truth_ugly-I_CN_ugly);
ErrU(3) = abs(truth_ugly-I_ugly(1));
ErrU(4) = abs(truth_ugly-I_ugly(2));
ErrU(5) = abs(truth_ugly-I_ugly(3));
ErrU    = ErrU./truth_ugly;

% Errors Exp
ErrE(1) = abs(truth_exp-q_exp);
ErrE(2) = abs(truth_exp-I_CN_exp);
ErrE(3) = abs(truth_exp-I_exp(1));
ErrE(4) = abs(truth_exp-I_exp(2));
ErrE(5) = abs(truth_exp-I_exp(3));
ErrE    = ErrE./truth_exp;

%% Plot

% Functions
figure(1)  
plot(x,f_ugly(x),'r-','Linewidth',2)
hold on
set(gca, 'FontName', 'Helvetica','FontSize',16)
xlabel('x')
ylabel('f(x)')
title('f(x) = x/2 + ((1/10 + x)sin(5x - 1))/(1 + x^2sin(x-0.5)^2)') 

figure(2)  
plot(x,f_exp(x),'r-','Linewidth',2)
hold on
set(gca, 'FontName', 'Helvetica','FontSize',16)
xlabel('x')
ylabel('f(x)')
title('f(x) = (4x + 4)exp(4x + 4)') 

% Errors
figure(3)
semilogy(1,(ErrU(3)),'b.','MarkerSize',30)
hold on
grid on
semilogy(2,(ErrU(4)),'r.','MarkerSize',30)
semilogy(3,(ErrU(5)),'g.','MarkerSize',30)
semilogy(4,(ErrU(1)),'m.','MarkerSize',30)
semilogy(5,(ErrU(2)),'k.','MarkerSize',30)
set(gca, 'FontName', 'Helvetica','FontSize',16)
xlabel('Method')
ylabel('Integral Error')
title('f(x) = x/2 + ((1/10 + x)sin(5x - 1))/(1 + x^2sin(x-0.5)^2)') 
legend('Gauss Two-point','Gauss Three-point','Gauss Four-point',...
    'MATLAB quad','Clenshaw-Curtis')

figure(4)
semilogy(1,(ErrE(3)),'b.','MarkerSize',30)
hold on
grid on
semilogy(2,(ErrE(4)),'r.','MarkerSize',30)
semilogy(3,(ErrE(5)),'g.','MarkerSize',30)
semilogy(4,(ErrE(1)),'m.','MarkerSize',30)
semilogy(5,(ErrE(2)),'k.','MarkerSize',30)
set(gca, 'FontName', 'Helvetica','FontSize',16)
xlabel('Method')
ylabel('Integral Error')
title('f(x) = (4x + 4)exp(4x + 4)') 
legend('Gauss Two-point','Gauss Three-point','Gauss Four-point',...
    'MATLAB quad','Clenshaw-Curtis')

% Function Evaluations
figure(5)
hold on
bar(1,2,'b')
bar(2,3,'r')
bar(3,4,'g')
bar(4,fevals_ugly,'m')
bar(5,M+1,'k')
set(gca, 'FontName', 'Helvetica','FontSize',16)
xlabel('Method')
ylabel('Function Evaluations')
title('f(x) = x/2 + ((1/10 + x)sin(5x - 1))/(1 + x^2sin(x-0.5)^2)') 
legend('Gauss Two-point','Gauss Three-point','Gauss Four-point',...
    'MATLAB quad','Clenshaw-Curtis','Location','NorthWest')

figure(6)
hold on
bar(1,2,'b')
bar(2,3,'r')
bar(3,4,'g')
bar(4,fevals_exp,'m')
bar(5,M+1,'k')
set(gca, 'FontName', 'Helvetica','FontSize',16)
xlabel('Method')
ylabel('Function Evaluations')
title('f(x) = (4x + 4)exp(4x + 4)') 
legend('Gauss Two-point','Gauss Three-point','Gauss Four-point',...
    'MATLAB quad','Clenshaw-Curtis','Location','NorthWest')

