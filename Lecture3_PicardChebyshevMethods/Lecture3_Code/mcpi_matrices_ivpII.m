% Texas A&M University - Department of Aerospace Engineering
% Author        : Robyn Woollands
% File name     : mcpi_matrices_ivpII.m
% Description   : Generates MCPI constant matrices for the second order IVP
% Date Written  : May 15, 2015
% Date Modified :
% References
%================================================================
function [Cxx,Cxv,Cxf,Ca,Cg,Cf] = mcpi_matrices_ivpII(N,M)
% INPUT         : N - order of Chebyshev Polynomial
%                 M - number of sample points
%               Notes
%                   M = N (interpolation)
%                   M > N (least squares)
% OUTPUT        : Cx, Ca matrices

% Chebyshev Polynomials
% R.H.S. VELOCITY
T_F = zeros(N-1,M+1);
for i = 0:N-2
    for j = 0:M
        T_F(i+1,j+1)    = cos(i*j*pi/M);
    end
end
% R.H.S. POSITION
T_beta = zeros(N,M+1);
for i = 0:N-1
    for j = 0:M
        T_beta(i+1,j+1)    = cos(i*j*pi/M);
    end
end

% L.H.S.
T_alpha = zeros(N+1,M+1);
for i = 0:N
    for j = 0:M
        T_alpha(i+1,j+1) = cos(i*j*pi/M);
    end
end
% Weight Matrix
W_tau                    = diag([0.5 ones(1,M-1) 0.5]);

% V Matrix
V                        = W_tau.*(2/M);

% R Matrix
% VELOCITY
R_beta                   = diag([1 1./(2.*[1:N-1])]);
% POSITION
if M == N; p = 1; end
if M > N;  p = 2; end
R_alpha                  = diag([1 1./(2.*[1:N-1]) 1/(p*N)]);

% S Matrix
% VELOCITY
Si                      = zeros(N,N-1);
Si(1,1)                 = 1;
Si(1,2)                 = -1/2;
for k = 3:N-1;
    Si(1,k)             = ((-1)^(k)) * ((1/(k - 2)) - (1/k)); % 1st row of S matrix
end
Sa                      = eye(N-1,N-1);
Sb                      = zeros(N-1,N-1);
Sb(1:end,3:end)         = -eye(N-1,N-3);
Si(2:end,1:end)         = Sa + Sb;
S_beta                  = Si;
clear Si
% POSITION
Si                      = zeros(N+1,N);
Si(1,1)                 = 1;
Si(1,2)                 = -1/2;
for k = 3:N;
    Si(1,k)             = ((-1)^(k)) * ((1/(k - 2)) - (1/k)); % 1st row of S matrix
end
Sa                      = eye(N,N);
Sb                      = zeros(N,N);
Sb(1:end,3:end)         = -eye(N,N-2);
Si(2:end,1:end)         = Sa + Sb;
S_alpha                 = Si;

% Output Matrices
Cf                      = T_F*V;
Ca                      = R_beta*S_beta*T_F*V;      % velocity
Cg                      = R_alpha*S_alpha;          % position

% CaII = R_alpha*S_alpha*R_beta*S_beta*T_F*V;
% CaI  = R_alpha*S_alpha;
if M == N; W_alpha = diag([0.5 ones(1,N-1) 0.5]);...
        W_beta = diag([0.5 ones(1,N-1)]); W_gamma = diag([0.5 ones(1,N-2)]);end
if M > N;  W_alpha = diag([0.5 ones(1,N)]); W_beta = diag([0.5 ones(1,N-1)]); end
Cxx                = T_alpha'*W_alpha;
Cxv                = T_beta'*W_beta;
Cxf                = T_F'*W_gamma;