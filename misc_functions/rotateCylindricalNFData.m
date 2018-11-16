function [data_nf_rotated] = rotateCylindricalNFData(data_nf,cylinder_axis)

% Create Results Table
p = length(data_nf.E);
data_nf_rotated = table(zeros(p,1),zeros(p,1),zeros(p,1),zeros(p,1),zeros(p,1),zeros(p,3),zeros(p,1));
data_nf_rotated.Properties.VariableNames = {'x','y','z','r','phi','E','Eabs'};

data_nf_rotated.x = data_nf.x;
data_nf_rotated.y = data_nf.y;
data_nf_rotated.z = data_nf.z;

switch cylinder_axis
    case 'x' 
        r = sqrt(data_nf.y.^2 + data_nf.z.^2);
        phi = atan2(data_nf.z,data_nf.y);
        for i =1:length(data_nf.E)
            R = [1 0 0;...
                 0 cos(phi(i)) sin(phi(i));...
                 0 -sin(phi(i)) cos(phi(i))];
            data_nf_rotated.E(i,:) = (R*data_nf.E(i,:)')';
            data_nf_rotated.r(i) = r(i);
            data_nf_rotated.phi(i) = phi(i);
            data_nf_rotated.Eabs(i) = data_nf.Eabs(i);
        end
        
    case 'y'
        r = sqrt(data_nf.x.^2 + data_nf.z.^2);
        phi = round(atan2(data_nf.x,data_nf.z),6);
        for i =1:length(data_nf.E)
            R = [cos(phi(i)) 0 sin(phi(i)) ;...
                 0 1 0;...
                -sin(phi(i)) 0 cos(phi(i))];
            data_nf_rotated.E(i,:) = (R*data_nf.E(i,:)')';
            data_nf_rotated.r(i) = r(i);
            data_nf_rotated.phi(i) = phi(i);
            data_nf_rotated.Eabs(i) = data_nf.Eabs(i);
        end

    case 'z' 
        r = sqrt(data_nf.x.^2 + data_nf.y.^2);
        phi = round(atan2(data_nf.y,data_nf.x),6);
        idx = find(phi < 0);
        phi(idx) = 2*pi + phi(idx);
        for i =1:length(data_nf.E)
            R = [cos(phi(i)) sin(phi(i)) 0;...
                 -sin(phi(i)) cos(phi(i)) 0;...
                 0 0 1];
            data_nf_rotated.E(i,:) = (R*data_nf.E(i,:)')';
            data_nf_rotated.r(i) = r(i);
            data_nf_rotated.phi(i) = phi(i);
            data_nf_rotated.Eabs(i) = data_nf.Eabs(i);
        end
end

