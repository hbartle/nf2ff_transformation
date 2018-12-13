function [xtheta,xphi] = sphericalVectorWaveFunctionFarField(s,m,n,A,theta,phi,k)
%SPHERICALVECTORWAVEFUNCTION
% Far-Field approximation of the spherical vector wave function
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
%       xtheta  Theta Component of Spherical wave function
%
%       xphi    Phi Component of Spherical wave function
%


% Helper Functions
delta = @(s,sigma) 0.5*(1+(-1)^(s+sigma));


%%%% Note %%%
% T and P provide values for each |m| to speed up the process
% This means the f output will be a matrix with m columns

% Theta Function
T = @(s,m,n,theta) delta(s,1)*1j*m./sin(theta) .* associatedLegendre(n,cos(theta))...
                   + delta(s,2).*associatedLegendreDerivative(n,cos(theta));
% Phi Function
P = @(m,phi) exp(1j*m.*repmat(phi,1,size(m,2)));




C = 1/(k*A) * (-1j)^(n+delta(s,1)) * exp(1j*k*A) *P(m,phi);
% Theta Component
xtheta = C .* T(s,m,n,theta);
% Phi Component
xphi = C .* (-1)^s .* T(s+1,m,n,theta);




end

