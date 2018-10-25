function [data_nf2ff] = nf2ffTransformation(data_nf,f,phi_range,theta_range)


% Wave Number
k = 2*pi*f/physconst('LightSpeed');

% Far-Field Point Radius (1m)
r = 1;


p = length(phi_range);
t = length(theta_range);

% Create Table
data_nf2ff = table(zeros(p*t,1),zeros(p*t,1),zeros(p*t,3),zeros(p*t,1));
data_nf2ff.Properties.VariableNames = {'theta','phi','E','Eabs'};

n = 1;
for phi=phi_range
    for theta=theta_range
        % Electric Field Contributions from each NF-measurement
        E_est = data_nf.E.*...
                exp(1j*k*(data_nf.z*cos(theta*pi/180)+...
                          data_nf.y*sin(theta*pi/180)*cos(phi*pi/180)+...
                          data_nf.x*sin(theta*pi/180)*sin(phi*pi/180)));
        % Sum Up Contributions
        E_est = exp(-1j*k*r)/r*sum(E_est);
        
        % Fill Table
        data_nf2ff.E(n,:) = E_est;
        data_nf2ff.Eabs(n) = sqrt(E_est*E_est');
        data_nf2ff.theta(n) = theta;
        data_nf2ff.phi(n) = phi;
        %disp(data_nf2ff.Eabs(n))
        n = n + 1;
    end
end
end


