% Robyn Woollands 2017
% Texas A&M University - Department of Aerospace Engineering
% File name     : clenshaw_curtis_ivpII.m
% Description   : Generates constant matrices for first order 
%                 Clenshaw & Curtis Quadrature
% Date Written  : March 27, 2017
% Date Modified : March 27, 2017
%
% Input:  N   -- Chebyshev polynomial order
%         M   -- Number of sample points
%
% Output: T2  -- Chebyshev Matrix [(M+1)x(N+1)]
%         P2  -- Picard Integration Operator [(N+1)xN]
%         T1  -- Chebyshev Matrix [(M+1)xN]
%         P1  -- Picard Integration Operator [Nx(N-1)]
%         Ta  -- Chebyshev Matrix [(M+1)x(N-1)]
%         A   -- Least Squares Operator [(N-1)x(M+1)]
%
% References: 
%================================================================

function [T2,P2,T1,P1,Ta,A] = clenshaw_curtis_ivpII(N,M)

% Least Squares Operator (A) [(N-1)x(M+1)]
s = -1;
[A,Ta] = lsq_chebyshev_fit(s,N-2,M);

% Compute "Velocity" Constants of Integration (i.e. evaluated T at tau = -1).
Lv = zeros(N,N); % [NxN]
for k = 0:N-1
    Lv(1,k+1) = cos(k*acos(-1)); % Const of Integration (CoI)
end

% Compute "Position" Constants of Integration (i.e. evaluated T at tau = -1).
Lp = zeros(N+1,N+1); % [NxN]
for k = 0:N
    Lp(1,k+1) = cos(k*acos(-1)); % Const of Integration (CoI)
end

% S Matrix for Velocity (Sv) [Nx(N-1)]
temp1              = diag([1 1./(2.*(1:N-1))]);
temp2              = eye(N-1,N-1);
temp3              = zeros(N-1,N-1);
temp3(1:end,3:end) = -eye(N-1,N-3);
temp4              = zeros(N,N-1);
temp4(2:end,1:end) = temp2 + temp3;
Sv                 = temp1*temp4;
Sv(1,1)            = 0.25;
Sv(2,1)            = 1;

% Picard Integration Operator (Acceleration to Velocity)
P1     = (eye(N,N) - Lv)*Sv;    % [Nx(N-1)]

clear temp1 temp2 temp3 temp4
% S Matrix for Position (Sp) [(N+1)xN]
temp1              = diag([1 1./(2.*(1:N))]);
temp2              = eye(N,N);
temp3              = zeros(N,N);
temp3(1:end,3:end) = -eye(N,N-2);
temp4              = zeros(N+1,N);
temp4(2:end,1:end) = temp2 + temp3;
Sp                 = temp1*temp4;
Sp(1,1)            = 0.25;
Sp(2,1)            = 1;

% Picard Integration Operator (Velocity to Position)
P2     = (eye(N+1,N+1) - Lp)*Sp;    % [(N+1)xN]

% Chebyshev Matrix (interpolate "velocity" coefficients)
[T1,~] = chebyshev(s,N-1,M,2);    % arg3 = 2 -> Trig Cheby Poly

% Chebyshev Matrix (interpolate "position" coefficients)
[T2,~] = chebyshev(s,N,M,2);      % arg3 = 2 -> Trig Cheby Poly
return
