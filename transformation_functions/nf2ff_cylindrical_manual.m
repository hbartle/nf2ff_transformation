function [data_nf2ff] = nf2ff_cylindrical_manual(data_nf,f,theta_range,window)

Ez = data_nf.E(:,3);
Ephi = data_nf.E(:,2);

% Wave Number
lambda = physconst('LightSpeed')/f;
k0 = 2*pi/lambda;

   
% Cylinder Radius
r0 = mean(data_nf.r);

% Number of Samples in z dimension
M = numel(unique(data_nf.z));
% Number of Samples in phi dimension
N = length(data_nf.phi)/M;

% Step size height dimension
delta_z = (max(data_nf.z)-min(data_nf.z))/(M-1);
% Step size angular dimension
delta_phi = 2*pi/N;

% Calculate k vector for each far field direction
n = 1;
for phi= unique(data_nf.phi)'
    for theta = theta_range
        k(n,:) = k0*[sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
       
        % Fill table with angle values
        data_nf2ff.theta(n) = theta;
        data_nf2ff.phi(n) = phi;
        n =n+1;
    end
end

% Hankel Function of second kind and its derivative
H = @(nu,Z) besselh(nu,2,Z);
dH = @(nu,Z) 0.5*(besselh(nu-1,2,Z) - besselh(nu+1,2,Z));



% Compute Helper Matrices
[n_grid, k_grid] = meshgrid(-N/2:N/2-1,k(:,2));
n_NM = reshape(n_grid,numel(n_grid),1);
k_NM = reshape(k_grid,numel(k_grid),1);
alpha_NM = sqrt(k0^2 - k_NM.^2);
% Hankel Coefficients
h = 1/(k0*r0)*n_NM.*k_NM .* H(n_NM,alpha_NM*r0);
dh = -dH(n_NM,alpha_NM*r0);
h2 = 1/k0 * alpha_NM .* H(n_NM,alpha_NM*r0);

A1 = exp(-1j*n_NM*data_nf.phi');
A2 = exp(1j*k_NM*data_nf.z');
A3 = 1j.^repmat(n_NM,1,N*M).*exp(1j*n_NM*repmat(data_nf.phi,M,1)');

b = 1./(alpha_NM/k0.*h).*(A2.*A1)*Ez;
a = 1./(dh).*(b.*n_NM.*k_NM/(k0*r0).*h - (A2.*A1)*Ephi);

% Far-Field Amplitude Factor
r = 100;
C = -2*k0*exp(-1j*k0*r)/r;
% Compute Far-Field


% Create Results Table
p = N;
t = length(theta_range);

data_nf2ff = table(zeros(p*t,1),zeros(p*t,1),zeros(p*t,1),zeros(p*t,1),zeros(p*t,1));
data_nf2ff.Properties.VariableNames = {'theta','phi','Etheta','Ephi','Eabs'};





end


