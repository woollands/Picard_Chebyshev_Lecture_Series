% Donghoon Kim & Robyn Woollands 2013
% Texas A&M University - Department of Aerospace Engineering
% File name     : FnG.m
% Description   : Computes Keplerian position and velocity at a given time 
% Date Written  : December 2, 2013
% Date Modified : March 26, 2017
%
% Input:  t0  -- Initial time
%         t   -- Current time
%         r0  -- Initial position vector
%         v0  -- Initial velocity vector
%         GM  -- Gravitational parameter
%
% Output: r    -- Position vector at time (t)
%         v    -- Velocity vector at time (t)
%         Ehat -- Eccentrical anomaly at time (t)
%
% References: 1. Schaub & Junkins, Analytical Mechanics of Space Systems (2010)
%================================================================

function [ r, v, Ehat ] = FnG(t0, t, r0, v0, GM)

    % Transform Row Vector to Column Vector
    if isrow(r0); r0 = r0'; end
    if isrow(v0); v0 = v0'; end

    % Find Mean Anomlay
    R0     = sqrt(r0'*r0);          % Initial position magnitude
    V0     = sqrt(v0'*v0);          % Initial velocity magnitude
    sigma0 = (r0'*v0)/sqrt(GM);     % Defined 
    A      = 2/R0 - V0^2/GM;        % Reciprocal of 1/a
    a      = 1/A;                   % Semimajor axis
    M      = sqrt(GM/a^3)*(t - t0); % Mean anomaly (rad)

    % Run Newton-Raphson Method 
    tol    = 1e-8;                  % Tolerance
    itr    = 0;                     % Initial iteration number
    MaxIt  = 10;                    % Maximum iteration number
    Ehat   = M;                     % Initial guess for eccentric anomaly (rad)
    dEhat  = 1;                     % Initial eccentric anomaly error (rad)
    while abs(dEhat) > tol;     
        err   = M - (Ehat - (1 - R0/a)*sin(Ehat) + sigma0/sqrt(a)*(1 - cos(Ehat)));
        derr  = - 1 + (1 - R0/a)*cos(Ehat) - sigma0/sqrt(a)*sin(Ehat);
        Ehat  = Ehat - err/derr;
        dEhat = - err/derr;
        itr   = itr + 1;
        if itr > MaxIt, break, end    
    end

    % Generate F & G Solution
    R      = a + (R0 - a)*cos(Ehat) + sqrt(a)*sigma0*sin(Ehat);
    F      = 1 - a/R0*(1 - cos(Ehat));
    G      = (t - t0) + sqrt(a^3/GM)*(sin(Ehat) - Ehat);
    Fdot   = - sqrt(GM*a)/(R*R0)*sin(Ehat);
    Gdot   = 1 - a/R*(1 - cos(Ehat));
    r      = F*r0 + G*v0;
    v      = Fdot*r0 + Gdot*v0;
    
return