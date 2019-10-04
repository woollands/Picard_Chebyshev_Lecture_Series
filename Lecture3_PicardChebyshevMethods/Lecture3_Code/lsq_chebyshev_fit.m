% Robyn Woollands 2013
% Texas A&M University - Department of Aerospace Engineering
% File name     : chebyshev.m
% Description   : Generates Chebyshev polynomials of the first kind
% Date Written  : June 1, 2013
% Date Modified : March 16, 2017
%
% Input:  s   -- sign on tau (-1 or 1)
%         N   -- Chebyshev polynomial order
%         M   -- Number of sample points
%
% Output: A   -- Least Squares Operator
%         T   -- Chebyshev matrix
%================================================================

function [A,T] = lsq_chebyshev_fit(s,N,M)

% Generate Chebyshev Polynomials
[T,~] = chebyshev(s,N,M,2);   % Default Trig Cheby

% Weight Matrix
W     = diag([0.5 ones(1,M-1) 0.5]);

% V Matrix
if M == N
    V     = diag([1/M (2/M).*ones(1,N-1) 1/M]);
elseif M > N
    V     = diag([1/M (2/M).*ones(1,N)]);
end

% Least Squares Operator
A = V*T'*W;

return