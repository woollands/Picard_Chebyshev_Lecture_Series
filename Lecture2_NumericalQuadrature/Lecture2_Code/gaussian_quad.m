% Robyn Woollands 2017
% Texas A&M University - Department of Aerospace Engineering
% File name     : gaussian_quad.m
% Description   : Computes the Gaussian Quadrature
% Date Written  : March 20, 2017
% Date Modified : March 20, 2017
%
% Input:  method -- Determines which quadrature method to use
%         f      -- Function to be integrated
%
% Output: I      -- Integral (area under the curve)
%================================================================

function I = gaussian_quad(method,f)

if method == 2
    % Two-point Formula
    pts = [-1/sqrt(3) 1/sqrt(3)];
    I   = sum(f(pts));
elseif method == 3
    % Three-point Formula
    pts = [-sqrt(3/5) 0 sqrt(3/5)];
    c1  = 5/9; 
    c2  = 8/9;
    I   = sum(c1*f(pts(1)) + c2*f(pts(2)) + c1*f(pts(3)));
elseif method == 4
    % Four-point Formula
    pts = [-0.861136 0.861136 -0.339981 0.339981];
    c1  = 0.34785;
    c2  = 0.652145;
    I   = sum(c1*(f(pts(1)) + f(pts(2))) + c2*(f(pts(3)) + f(pts(4))));
end

return