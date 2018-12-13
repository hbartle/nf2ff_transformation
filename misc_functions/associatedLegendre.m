function [P] = associatedLegendre(n,x)
%ASSOCIATEDLEGENDRE
% Modified associated legendre function to duplicate negative m values
%
%   Input Arguments:
%       
%           n       Degree
%           x       Argument
%
%   Output Arguments
%
%           P       Polynomial Value

p = legendre(n,x);
P = [flipud(p(2:end,:));p]' ;

end

