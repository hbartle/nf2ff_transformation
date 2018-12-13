function [fr,ftheta,fphi, Y] = sphericalVectorWaveFunction(s,m,n,A,theta,phi,k)
%SPHERICALVECTORWAVEFUNCTION
%
% Input Arguments:
%
%       s       Parity
%       m,n     Iterators over 
%       A       Radius
%       theta   Theta angle
%       phi     Phi angle
%       k       Wavenumber at given frequency
%
% Output Arguments:
%
%       fr      Radial Component of Spherical wave function
%
%       ftheta  Theta Component of Spherical wave function
%
%       fphi    Phi Component of Spherical wave function
%
%       Y       Factor that shows up when solving for the spherical wave
%               expansion coefficients
%

% Helper Functions
z =@(n,r) sqrt(pi/(2*r)).* besselh(n+0.5,r);
dz = @(n,r) z(n-1,r) - (n+1)/r * z(n,r); % CHECK IF THIS IS CORRECT!!
delta = @(s,sigma) 0.5*(1+(-1)^(s+sigma));


%%%% Note %%%
% R,T and P provide values for each |m| to speed up the process
% This means the f output will be a matrix with m columns

% Radial Function
R = @(s,n,r) delta(s,1).*z(n,k*r) + delta(s,2).*1/(k*r).* dz(n,k*r);
% Theta Function
T = @(s,m,n,theta) delta(s,1)*1j*m./sin(theta) .* associatedLegendre(n,cos(theta))...
                   + delta(s,2) .* associatedLegendreDerivative(n,cos(theta)); % Add derivative of Associated Legendre Function;

% Phi Function
P = @(m,phi) exp(1j*m.*repmat(phi,1,size(m,2)));



% Radial Component
fr = delta(s,2)*n*(n+1)/(k*A) *z(n,k*A) * associatedLegendre(n,cos(theta)).*P(m,theta);
% Theta Component
ftheta = repmat(R(s,n,A),1,size(m,2)).*T(s,m,n,theta).*P(m,phi);
% Phi Component
fphi = (-1)^s.*repmat(R(s,n,A),1,size(m,2)).*T(s+1,m,n,theta).*P(m,theta);



% Additional Factor that shows up when solving for the spherical wave
% expansion coefficients.
Y = 4*pi/(2*n + 1) *factorial(n+abs(m(1,:)))./factorial(n-abs(m(1,:))) *n*(n+1)*R(s,n,A).^2;

end

