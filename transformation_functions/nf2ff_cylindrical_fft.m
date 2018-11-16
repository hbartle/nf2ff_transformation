function [data_nf2ff] = nf2ff_cylindrical_fft(data_nf,f,theta_range,window)

% Wave Number
lambda = physconst('LightSpeed')/f;
k0 = 2*pi/lambda;

% Hankel Function of second kind and its derivative
H = @(nu,Z) besselh(nu,2,Z);
dH = @(nu,Z) 0.5*(besselh(nu-1,2,Z) - besselh(nu+1,2,Z));
   
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

n = -N/2:(N/2-1);
m = -M/2:(M/2-1);
delta_kz =2*pi/(M*delta_z); 
kz = m*delta_kz;

[k_grid, n_grid] = meshgrid(kz,n);

% Hankel Coefficients
h = 1/(k0*r0)*n_grid.*k_grid .* H(n_grid,sqrt(k0^2 - k_grid.^2)*r0);
dh = -dH(n_grid,sqrt(k0^2 - k_grid.^2)*r0);
h2 = 1/k0 * sqrt(k0^2 - k_grid.^2) .* H(n_grid,sqrt(k0^2 - k_grid.^2)*r0);

% Reshape Efield components
Ephi = reshape(data_nf.E(:,2),N,M);
Ez = reshape(data_nf.E(:,3),N,M);

% Sampling Window Function
[Lphi,Lz] = size(Ephi);
switch window
    case 'none'
        w_phi = ones(Lphi,1);
        w_z = ones(Lz,1);
    case 'tukey'
        w_phi = tukeywin(Lphi);
        w_z = tukeywin(Lz);
    case 'hamming'
        w_phi = hamming(Lphi);
        w_z = hamming(Lz);
    otherwise
        error('No valid window function selected')
        
end

W = w_phi*w_z'/(Lphi*Lz);

% Spectral analysis
% b = 1./(h2.*delta_kz).*fftshift(fft(ifftshift(ifft(Ez,[],1)),[],2));
% a = 1./(dh.*delta_kz).*(b.*h.*delta_kz - fftshift(fft(ifftshift(ifft(Ephi,[],1)),[],2)));

b = 1./(h2.*delta_kz).*fftshift(fft(ifft(Ez.*W,[],1),[],2));
a = 1./(dh.*delta_kz).*(b.*h.*delta_kz - fftshift(fft(ifft(Ephi.*W,[],1),[],2)));


% Spherical Far-Field Wavenumber Vector
[theta_grid,n_grid_spherical]=meshgrid(theta_range,n);
kz_grid_spherical = k0*cos(theta_grid);

a_ff=interp2(n_grid',k_grid',a',n_grid_spherical',kz_grid_spherical','spline')';
b_ff=interp2(n_grid',k_grid',b',n_grid_spherical',kz_grid_spherical','spline')';


% Far-Field Amplitude Factor
r = 100;
C = -2*k0*exp(-1j*k0*r)/r;
% Compute Far-Field
% Etheta = 1j*C*sqrt(k0^2-kz_grid_spherical.^2).*fftshift(fft(1j.^n_grid_spherical.*b_ff));
% Ephi = C*sqrt(k0^2-kz_grid_spherical.^2).*fftshift(fft(1j.^n_grid_spherical.*a_ff));
% Eabs = abs(sqrt(Etheta.^2 + Ephi.^2));


Etheta = 1j*C*sqrt(k0^2-kz_grid_spherical.^2).*ifft(1j.^n_grid_spherical.*b_ff);
Ephi = C*sqrt(k0^2-kz_grid_spherical.^2).*ifft(1j.^n_grid_spherical.*a_ff);
Eabs = abs(sqrt(Etheta.^2 + Ephi.^2));

% Create Results Table
p = N;
t = length(theta_range);

data_nf2ff = table(zeros(p*t,1),zeros(p*t,1),zeros(p*t,1),zeros(p*t,1),zeros(p*t,1));
data_nf2ff.Properties.VariableNames = {'theta','phi','Etheta','Ephi','Eabs'};

s = numel(theta_grid);
data_nf2ff.theta = reshape(theta_grid,s,1);
data_nf2ff.phi = reshape(delta_phi*(n_grid_spherical+abs(min(n_grid_spherical))),s,1);
data_nf2ff.Etheta = reshape(Etheta,s,1);
data_nf2ff.Ephi = reshape(Ephi,s,1);
data_nf2ff.Eabs = reshape(Eabs,s,1);


end


