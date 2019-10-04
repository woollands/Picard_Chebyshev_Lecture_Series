% Robyn Woollands 2015
% Texas A&M University - Department of Aerospace Engineering
% File name     : chebyshev.m
% Description   : Generates Chebyshev polynomials of the first kind
% Date Written  : June 1, 2013
% Date Modified : March 16, 2017
%
% Input:  s   -- sign on tau (-1 or 1)
%         N   -- Chebyshev polynomial order
%         M   -- Number of sample points
%         arg -- Recursive OR Trigonometric formulation
%
% Output: T   -- Chebyshev polynomials
%         tau -- Cosine Sample Points
%================================================================

function [T,tau] = chebyshev(s,N,M,arg)

% Cosine Sample Points
tau = s.*cos([0:M].*pi/M);

if arg == 1;
    
    % Chebyshev Polynomials (Recursive Formulation)
    T(:,1) = ones(M+1,1);
    T(:,2) = tau';
    for k = 2:N;
        T(:,k+1) = 2.*tau'.*T(:,k) - T(:,k-1);
    end
    
elseif arg == 2
    
    % Chebyshev Polynomials (Trigonometric Formulation)
    for k = 0:N
        T(:,k+1) = cos(k*acos(tau'));
    end
    
end

return