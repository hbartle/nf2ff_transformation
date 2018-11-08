function [data_nf2ff] = nf2ffTransformation_manual(data_nf,f,phi_range,theta_range,window)


% Create Results Table
p = length(phi_range);
t = length(theta_range);

data_nf2ff = table(zeros(p*t,1),zeros(p*t,1),zeros(p*t,3),zeros(p*t,1));
data_nf2ff.Properties.VariableNames = {'theta','phi','E','Eabs'};
    

% Wave Number
lambda = physconst('LightSpeed')/f;
k_const = 2*pi/lambda;
       
% Far-Field Point Radius [m]
%r = 2*pi;
r = 1;
% Calculate k vector for each far field direction
n = 1;
for phi= phi_range
    for theta = theta_range
        k(n,:) = k_const*[sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
        
        % Fill table with angle values
        data_nf2ff.theta(n) = theta;
        data_nf2ff.phi(n) = phi;
        
        n =n+1;
    end
end


% Get step size in x,y direction
delta_x = unique(data_nf.x);
delta_x = abs(delta_x(1)-delta_x(2));
delta_y = unique(data_nf.y);
delta_y = abs(delta_y(1)-delta_y(2));

% Sampling Window Function
L = length(unique(data_nf.x));

switch window
    case 'none'
        h = ones(L);
    case 'tukey'
        h = tukeywin(L);
    case 'hamming'
        h = hamming(L);
    otherwise
        disp('No valid window function selected')
end
        
H = reshape(h*h'/L,length(data_nf.Eabs),1);
Ex = data_nf.E(:,1).*H;
Ey = data_nf.E(:,2).*H;


% Create A matrix
% A = delta_x*delta_y * exp(1j*k*[data_nf.x,data_nf.y,data_nf.z]');
A = delta_x*delta_y * exp(1j*k(:,1:2)*[data_nf.x,data_nf.y]');

% Calculate Plane Wave Modal Coefficients
Fx = A*Ex;
Fy = A*Ey;
% Fz = A*data_nf.E(:,3);
Fz = -(k(:,1).*Fx + k(:,2).*Fy)./k(:,3);

% Calculate far-field E
data_nf2ff.E = 1j * exp(-1j*k_const*r)/(lambda*r)*1/k_const* repmat(k(:,3),1,3) .* [Fx Fy Fz];


data_nf2ff.Eabs = vecnorm(abs(data_nf2ff.E),2,2);


end


