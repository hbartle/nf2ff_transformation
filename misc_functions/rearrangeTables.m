function [data_nf] = rearrangeTables(data_nf)

Ex = data_nf.ExReal+ 1j*data_nf.ExImg;
Ey = data_nf.EyReal+ 1j*data_nf.EyImg;
Ez = data_nf.EzReal+ 1j*data_nf.EzImg;
Eabs = data_nf.EabsReal + 1j*data_nf.EabsImg;
%data_nf = table(data_nf.X,data_nf.Y,data_nf.Z,data_nf.R,data_nf.Theta,data_nf.Phi,[Ex,Ey,Ez],Eabs);
%data_nf.Properties.VariableNames = {'x' 'y' 'z' 'r' 'theta' 'phi' 'E' 'Eabs'};
data_nf = table(data_nf.X,data_nf.Y,data_nf.Z,[Ex,Ey,Ez],Eabs);
data_nf.Properties.VariableNames = {'x' 'y' 'z' 'E' 'Eabs'};

data_nf.x = data_nf.x/1000;
data_nf.y = data_nf.y/1000;
data_nf.z = data_nf.z/1000;


end

