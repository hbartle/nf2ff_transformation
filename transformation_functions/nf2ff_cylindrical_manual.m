function [data_nf2ff] = nf2ff_cylindrical_manual(data_nf,f,theta_range,phi_range,window)

% Hankel Function of second kind and its derivative
H = @(nu,Z) besselh(nu,2,Z);
dH = @(nu,Z) 0.5*(besselh(nu-1,2,Z) - besselh(nu+1,2,Z));

% Wave Number
lambda = physconst('LightSpeed')/f;
k0 = 2*pi/lambda;

Ez = data_nf.E(:,3);
Ephi = data_nf.E(:,2);

% Cylinder Radius
r0 = mean(data_nf.r)-0.15;
% Number of Samples in z dimension
M = numel(unique(data_nf.z));
% Number of Samples in phi dimension
N = length(data_nf.phi)/M;

% Step size height dimension
delta_z = (max(data_nf.z)-min(data_nf.z))/(M-1);
% Step size angular dimension
delta_phi = 2*pi/N;


% Calculate kz
kz = k0*cos(theta_range);
delta_kz = diff(kz);
delta_kz = [2*delta_kz(1)-delta_kz(2) delta_kz];
delta_kz =repelem(2*pi/(M*delta_z),length(kz)); 

% Sampling Window Function
[Lphi,Lz] = size(Ephi);
switch window
    case 'none'
        w_phi = ones(Lphi,1);
        w_z = ones(Lz,1);
    case 'tukey'
        w_phi = tukeywin(Lphi,0.3);
        w_z = tukeywin(Lz,0.3);
    case 'hamming'
        w_phi = hamming(Lphi);
        w_z = hamming(Lz);
    otherwise
        error('No valid window function selected')
        
end

W = w_phi*w_z';

% Calculate spectral coefficients
n = -N/2:N/2-1;
for n_iter = 1:N
    for m=1:length(kz)
        alpha =sqrt(k0^2 - kz(m)); 
        b(n_iter,m) = 1/(alpha/k0*H(n(n_iter),alpha*r0)*delta_kz(m)) *...
                       (Ez.*W)'*exp(-1j*n(n_iter)*data_nf.phi + 1j*kz(m)*data_nf.z);
        a(n_iter,m) = 1/(alpha*dH(n(n_iter),alpha*r0)*delta_kz(m))...
                      *(b(n_iter,m)*n(n_iter)*kz(m)/(k0*r0)*H(n(n_iter),alpha*r0)*delta_kz(m)...
                      -(Ephi.*W)'*exp(-1j*n(n_iter)*data_nf.phi + 1j*kz(m)*data_nf.z));
    end
end


% Far-Field Amplitude Factor
r = 100;
C = -2*k0*exp(-1j*k0*r)/r;
% Compute Far-Field
for n_iter=1:length(phi_range)
    for m=1:length(kz)
        alpha =sqrt(k0^2 - kz(m)); 
        Etheta_ff(n_iter,m) = 1j*alpha*C* (1j.^n.*b(:,m)')*exp(1j*n*phi_range(n_iter))';
        Ephi_ff(n_iter,m) = alpha*C* (1j.^n.*a(:,m)')*exp(1j*n*phi_range(n_iter))';
       
        Eabs(n_iter,m) = abs(sqrt(Etheta_ff(n_iter,m)^2 + Ephi_ff(n_iter,m)^2));
    end
end

% Rearrange E-fields to have proper alignment in phi dimension
Etheta_ff = [Etheta_ff(end/2+1:end,:);Etheta_ff(1:end/2,:)];
Ephi_ff = [Ephi_ff(end/2+1:end,:);Ephi_ff(1:end/2,:)];
Eabs = [Eabs(end/2+1:end,:);Eabs(1:end/2,:)];


% Create Results Table
s = numel(Etheta_ff);
data_nf2ff = table(zeros(s,1),zeros(s,1),zeros(s,1),zeros(s,1),zeros(s,1));
data_nf2ff.Properties.VariableNames = {'theta','phi','Etheta','Ephi','Eabs'};

[t,p] = meshgrid(acos(kz/k0),phi_range);
data_nf2ff.theta = reshape(t,s,1);
data_nf2ff.phi = reshape(p,s,1);
data_nf2ff.Etheta = reshape(Etheta_ff,s,1);
data_nf2ff.Ephi = reshape(Ephi_ff,s,1);
data_nf2ff.Eabs = reshape(Eabs,s,1);



end


