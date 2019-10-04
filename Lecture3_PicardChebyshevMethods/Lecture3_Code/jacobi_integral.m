%================={ Copyright (c) Ahmad Bani-Younes 2013 }=================
%                Texas A&M University - AEROSPACE department
% 
% File name     : Jacobi_Integral.m
% Subject       : Computes the Hamiltonian.
% Description   : To validate orbit propagators accuracy.
%    
% Sub-files     : egm2008GPsphericalharmonic.m
%
% Compiler      : MATLAB 7.11.0 (R2010b)
% Date Modified : 03/27/2017
%==========================================================================

function dEnrgy = jacobi_integral(soln)

global GM Re omega Deg C S

load('aeroegm2008.mat') % [GM, Re, degree, C, S]

r_B    = soln(:,1:3);       % position in the rotating frame
v_B    = soln(:,4:6);       % velocity in the rotating frame
term   = 0.5*omega*omega.*(r_B(:,1).^2 + r_B(:,2).^2);
KE     = 0.5*( v_B(:,1).^2 + v_B(:,2).^2 + v_B(:,3).^2 );
V      = -egm2008GPsphericalharmonic(r_B, Deg) ;
Enrgy  = V + KE - term;
Enrgyo = Enrgy(1); 
dEnrgy = abs((Enrgy - Enrgyo)/Enrgyo);

return