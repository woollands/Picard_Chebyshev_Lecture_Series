function acc = ecef2eci(t,aB,omega)
% ecef2eci.m
% AUTHOR:       Robyn Woollands
% DATE:         June 15, 2016
% DESCRIPTION:  Function to convert acceleration from body to inertial
% frame

% Inputs:   t  -- time (s)
%          axB -- ECEF Position (km)
%        omega -- Earth rotational velocity (rad/s)
%
% Outputs: acc -- ECEF Acceleration (km/s^2)

th           = t*omega;
cos_th       = cos(th);
sin_th       = sin(th);

% Convert to Acceleration to Inertial Frame
acc(:,1)     = cos_th.*aB(:,1) - sin_th.*aB(:,2);
acc(:,2)     = sin_th.*aB(:,1) + cos_th.*aB(:,2);
acc(:,3)     = aB(:,3);

return