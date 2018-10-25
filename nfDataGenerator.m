%
% Creates a Sine-Shaped E field to simulate a waveguide opening
%
close all
% Get the X values for the setup
xValues = data_nf{1}.x;

% Get Spatial Frequency
k = 2*pi/(2*(max(xValues)-min(xValues)));

% E-Field  
E_norm = sin(k*xValues + pi/2 );

data_nf{1}.E = [zeros(length(xValues),1), E_norm, zeros(length(xValues),1)];

data_nf{1}.Eabs = vecnorm(data_nf{1}.E,2,2);