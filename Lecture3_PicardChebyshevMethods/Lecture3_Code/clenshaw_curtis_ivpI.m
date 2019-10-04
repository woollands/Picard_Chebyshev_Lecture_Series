% Robyn Woollands 2017
% Texas A&M University - Department of Aerospace Engineering
% File name     : clenshaw_curtis_ivpI.m
% Description   : Generates constant matrices for first order 
%                 Clenshaw & Curtis Quadrature
% Date Written  : March 26, 2017
% Date Modified : March 26, 2017
%
% Input:  N   -- Chebyshev polynomial order
%         M   -- Number of sample points
%
% Output: T1  -- Chebyshev Matrix [(M+1)x(N+1)]
%         P1  -- Picard Integration Operator [(N+1)xN]
%         Ta  -- Chebyshev Matrix [(M+1)xN]
%         A   -- Least Squares Operator [Nx(M+1)]
%
% References: 
%================================================================

function [T1,P1,Ta,A] = clenshaw_curtis_ivpI(N,M)

% Least Squares Operator (A)
s      = -1;    % Sign on tau
[A,Ta] = lsq_chebyshev_fit(s,N-1,M);

% Compute Constants of Integration (i.e. evaluated T at tau = -1).
L = zeros(N+1,N+1);
for k = 0:N
    L(1,k+1) = cos(k*acos(-1)); % Const of Integration (CoI)
end

% S Matrix
temp1              = diag([1 1./(2.*(1:N))]);
temp2              = eye(N,N);
temp3              = zeros(N,N);
temp3(1:end,3:end) = -eye(N,N-2);
temp4              = zeros(N+1,N);
temp4(2:end,1:end) = temp2 + temp3;
S                  = temp1*temp4;
S(1,1)             = 0.25;
S(2,1)             = 1;

% Picard Integration Operator
P1     = (eye(N+1,N+1) - L)*S;

% Chebyshev Matrix (interpolate solution coefficients)
[T1,~] = chebyshev(s,N,M,2);    % arg3 = 2 -> Trig Cheby Poly

return
