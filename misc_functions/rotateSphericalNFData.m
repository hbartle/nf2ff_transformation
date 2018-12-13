function [data_nf_rotated] = rotateSphericalNFData(data_nf)

% Create Results Table
p = length(data_nf.E);
data_nf_rotated = table(zeros(p,1),zeros(p,1),zeros(p,1),zeros(p,1),zeros(p,1),zeros(p,1),zeros(p,3),zeros(p,1));
data_nf_rotated.Properties.VariableNames = {'x','y','z','r','theta','phi','E','Eabs'};

data_nf_rotated.x = data_nf.x;
data_nf_rotated.y = data_nf.y;
data_nf_rotated.z = data_nf.z;

r = sqrt(data_nf.x.^2 + data_nf.y.^2 + data_nf.z.^2);
theta = acos(data_nf.z./r);
phi = atan2(data_nf.y,data_nf.x);
for i =1:length(data_nf.E)
    R = [sin(theta(i))*cos(phi(i)) sin(theta(i))*sin(phi(i)) cos(theta(i));...
         cos(theta(i))*cos(phi(i)) cos(theta(i))*sin(phi(i)) -sin(theta(i));...
         -sin(phi(i)) cos(phi(i)) 0];
     
    %E(1)-> Radial,E(2)->Theta,E(3)->Phi
    data_nf_rotated.E(i,:) = (R*data_nf.E(i,:)')';
    data_nf_rotated.r(i) = r(i);
    data_nf_rotated.theta(i) = theta(i);
    if phi(i) <0
        data_nf_rotated.phi(i) = pi - phi(i);
    else
        data_nf_rotated.phi(i) = phi(i);
    end
    data_nf_rotated.Eabs(i) = data_nf.Eabs(i);
end
  