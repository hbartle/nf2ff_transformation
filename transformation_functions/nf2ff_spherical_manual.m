function [data_nf2ff] = nf2ff_spherical_manual(data_nf,f,theta_range,phi_range)

% Wave Number
lambda = physconst('LightSpeed')/f;
k0 = 2*pi/lambda;
% Antenna minimal sphere radius
r0 = 2/3*mean(data_nf.r);
% Radius of Measurement Sphere
A = mean(data_nf.r);

theta = data_nf.theta;
phi = data_nf.phi;
Etheta = data_nf.E(:,2);
Ephi = data_nf.E(:,3);

% Step Sizes
delta_theta = diff(unique(round(data_nf.theta,2)));
delta_theta = delta_theta(1);
delta_phi = diff(unique(round(data_nf.phi,2)));
delta_phi = delta_phi(1);

% Theta/Phi values of Far-field points to calculate
[theta_ff,phi_ff] = meshgrid(theta_range,phi_range);
theta_ff = reshape(theta_ff,numel(theta_ff),1);
phi_ff = reshape(phi_ff,numel(phi_ff),1);
% Far-Field radius
r_ff = 10;


% Number of spherical wave modes
N = round(k0*r0 + 5);

% Preallocate space for increased speed
Etheta_ff = zeros(size(theta_ff));
Ephi_ff = zeros(size(phi_ff));

idx = 1;
for s=1:2
    for n = 1:N
        M = repmat((-n:n),length(theta),1);
        % Compute Spherical Expansion Coefficients
        [~,ftheta,fphi,Y] = sphericalVectorWaveFunction(s,M,n,A,theta,phi,k0);
        ftheta_t = conj(ftheta).*sin(theta);
        fphi_t = conj(fphi).*sin(theta);
        q = (1./Y .*(Etheta'*ftheta_t + Ephi'*fphi_t))'*delta_theta*delta_phi;
        
        % Compute Far-Field 
        M = repmat((-n:n),length(theta_ff),1);
        [xtheta,xphi] = sphericalVectorWaveFunctionFarField(s,M,n,r_ff,theta_ff,phi_ff,k0);

        Etheta_ff = Etheta_ff + sum(repmat(q,1,length(theta_ff))'.*xtheta,2);
        Ephi_ff = Ephi_ff + sum(repmat(q,1,length(theta_ff))'.*xphi,2);

        idx = idx+1
    end
end


% Create Results Table
s = numel(Etheta_ff);
data_nf2ff = table(zeros(s,1),zeros(s,1),zeros(s,1),zeros(s,1),zeros(s,1));
data_nf2ff.Properties.VariableNames = {'theta','phi','Etheta','Ephi','Eabs'};


data_nf2ff.theta = reshape(t,s,1);
data_nf2ff.phi = reshape(p,s,1);
data_nf2ff.Etheta = reshape(Etheta_ff,s,1);
data_nf2ff.Ephi = reshape(Ephi_ff,s,1);
data_nf2ff.Eabs = reshape(Eabs,s,1);



end


