function [xB,vB] = eci2ecef(t,X,V,omega)
% eci2ecef.m
% AUTHOR:       Robyn Woollands
% DATE:         June 15, 2016
% DESCRIPTION:  Function to convert states to body frame

% Inputs:  t  -- time (s)
%          X  -- ECI Position (km)
%          V  -- ECI Velocity (km/s)
%      omega  -- Earth rotational velocity (rad/s)
%
% Outputs: xB -- ECEF Position (km)
%          vB -- ECEF Velocity (km/s)

th              = t*omega;
cos_th          = cos(th);
sin_th          = sin(th);

% Convert to Position to Rotating Frame
xB(:,1)   =  cos_th.*X(:,1) + sin_th.*X(:,2);
xB(:,2)   = -sin_th.*X(:,1) + cos_th.*X(:,2);
xB(:,3)   = X(:,3);
% Convert Velocity to the Rotating Frame
if sum(sum(V)) ~= 0;
    vB(:,1)   = V(:,1) + omega*X(:,2);
    vB(:,2)   = V(:,2) - omega*X(:,1);
    vB(:,3)   = V(:,3);
elseif sum(sum(V)) == 0;
    vB = 0;
end

return