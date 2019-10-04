% Robyn Woollands 2013
% Texas A&M University - Department of Aerospace Engineering
% File name     : chebyshev.m
% Description   : Generates Chebyshev polynomials of the first kind
% Date Written  : June 1, 2013
% Date Modified : March 16, 2017
%
% Input:  N   -- Chebyshev polynomial order
%         M   -- Number of sample points
%
% Output: A   -- Least Squares Operator [(N+1)x(M+1)]
%         T   -- Chebyshev matrix [(M+1)x(N+1)]
%================================================================

function [A,T] = lsq_chebyshev_fit(N,M)

% Generate Chebyshev Polynomials [(M+1)x(N+1)]
[T,~] = chebyshev(N,M,2);   % Default Trig Cheby

% Weight Matrix [(M+1)x(M+1)]
W     = diag([0.5 ones(1,M-1) 0.5]);

% V Matrix [(N+1)x(N+1)]
if M == N
    V     = diag([1/M (2/M).*ones(1,N-1) 1/M]);
elseif M > N
    V     = diag([1/M (2/M).*ones(1,N)]);
end

% Least Squares Operator [(N+1)x(M+1)]
A = V*T'*W;

return