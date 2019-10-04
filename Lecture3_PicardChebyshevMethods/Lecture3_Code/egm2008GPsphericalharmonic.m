function GravityPotential = egm2008GPsphericalharmonic( p, maxdeg )
%  Created by Ahmad 12/20/2012 'EGM2008' to calculate the gravity potential;
global GM  Re C S
% Compute geocentric radius
r = sqrt( sum( p.^2, 2 ));
% Compute geocentric latitude
phic = asin( p(:,3)./ r );
% Compute lambda                                                           
lambda   = atan2( p(:,2), p(:,1) );
smlambda = zeros( size(p,1), maxdeg );
cmlambda = zeros( size(p,1), maxdeg );

slambda = sin(lambda);
clambda = cos(lambda);
smlambda(:,1) = 0;
cmlambda(:,1) = 1;
smlambda(:,2) = slambda;
cmlambda(:,2) = clambda;

for m=3:maxdeg+1
    smlambda(:,m) = 2.0.*clambda.*smlambda(:, m-1) - smlambda(:, m-2);
    cmlambda(:,m) = 2.0.*clambda.*cmlambda(:, m-1) - cmlambda(:, m-2);
end

% Compute normalized associated legendre polynomials
[P scaleFactor] = loc_gravLegendre( phic, maxdeg );

% Compute gravity in ECEF coordinates
GravityPotential = loc_gravityPCPF( p, maxdeg, P, C( 1:maxdeg+1, 1:maxdeg+1 ), ...
                                  S( 1:maxdeg+1, 1:maxdeg+1 ), smlambda, ...
                                  cmlambda, GM, Re, r, scaleFactor );
end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
function [P scaleFactor] = loc_gravLegendre( phi, maxdeg )
% loc_GRAVLEGENDRE internal function computing normalized associated 
% legendre polynomials, P, via recursion relations for spherical harmonic
% gravity 

P = zeros(maxdeg+2, maxdeg+2, length(phi));
scaleFactor = zeros(maxdeg+3, maxdeg+3, length(phi));

cphi = cos(pi/2-phi);
sphi = sin(pi/2-phi);

% force numerically zero values to be exactly zero
cphi(abs(cphi)<=eps) = 0;
sphi(abs(sphi)<=eps) = 0;
 
% Seeds for recursion formula
P(1,1,:) = 1;            % n = 0, m = 0;
scaleFactor(1,1,:) = 0;
P(2,1,:) = sqrt(3)*cphi; % n = 1, m = 0;
scaleFactor(2,1,:) = 1;
P(2,2,:) = sqrt(3)*sphi; % n = 1, m = 1;
scaleFactor(2,2,:) = 0;

for n = 2:maxdeg+2
    k = n + 1;
    for m = 0:n
        p = m + 1;
        % Compute normalized associated legendre polynomials, P, via recursion relations 
        % Scale Factor needed for normalization of dUdphi partial derivative
                
        if (n == m)           
            P(k,k,:) = sqrt(2*n+1)/sqrt(2*n)*sphi.*reshape(P(k-1,k-1,:),size(phi));
            scaleFactor(k,k,:) = 0;
        elseif (m == 0)
            P(k,p,:) = (sqrt(2*n+1)/n)*(sqrt(2*n-1)*cphi.*reshape(P(k-1,p,:),size(phi)) - (n-1)/sqrt(2*n-3)*reshape(P(k-2,p,:),size(phi)));
            scaleFactor(k,p,:) = sqrt( (n+1)*(n)/2);
        else
            P(k,p,:) = sqrt(2*n+1)/(sqrt(n+m)*sqrt(n-m))*(sqrt(2*n-1)*cphi.*reshape(P(k-1,p,:),size(phi)) - sqrt(n+m-1)*sqrt(n-m-1)/sqrt(2*n-3)*reshape(P(k-2,p,:),size(phi)));
            scaleFactor(k,p,:) = sqrt( (n+m+1)*(n-m));
        end
    end
end
end

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
function GPot = loc_gravityPCPF(p,maxdeg,P,C,S,smlambda,cmlambda,GM,Re,r,scaleFactor)
% loc_GRAVITYPCPF internal function computing gravity in planet-centered
% planet-fixed (PCEF) coordinates using PCPF position, desired
% degree/order, normalized associated legendre polynomials, normalized
% spherical harmonic coefficients, trigonometric functions of geocentric
% latitude and longitude, planetary constants, and radius to center of
% planet. Units are MKS.

rRatio   = Re./r;
rRatio_n = rRatio;

% initialize summation of gravity in radial coordinates
USumN         = 1;

% summation of gravity in radial coordinates
for n = 2:maxdeg
    k = n+1;
    rRatio_n      = rRatio_n.*rRatio;
    USumM         = 0;
    for m = 0:n
        j = m+1;
        USumM      = USumM + reshape(P(k,j,:),size(r)).*(C(k,j).*cmlambda(:,j) + S(k,j).*smlambda(:,j)); 
    end
    USumN         = USumN         + USumM.*rRatio_n;
end

% gravity in spherical coordinates
GPot       =  GM./r.*USumN;
end
