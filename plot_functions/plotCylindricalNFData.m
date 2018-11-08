function [] = plotCylindricalNFData(data_nf)

figure('name','Near-Field Magnitude','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1]);

h = unique(data_nf.y);
phi_half = unique(round(acos(data_nf.z/max(data_nf.z)),4));
phi = [phi_half(1:end-1); pi+phi_half(1:end-1)];
numelHeight= length(h);
numelPhi= length(phi);


[h_grid, phi_grid] = meshgrid(h,phi); 

Eabs =  reshape(data_nf.Eabs,numelPhi,numelHeight);
surf(phi_grid*180/pi,h_grid, Eabs )

title({'Cylindrical Scan - Magnitude';['Number of Probes in Linear Dimension: ' num2str(numelHeight)];...
       ['Number of Probes in Circular Dimension: ' num2str(numelPhi)]})
ylabel('Linear Dimension [mm]')
xlabel('Circular Dimension [°]')
zlabel('E-Field [V/m]')
% figure('name','Near-Field Phase','numbertitle','off',...
%         'units','normalized','outerposition',[0 0 1 1]);
% 
% numelX= sqrt(length(data_nf.x));
% numelY= sqrt(length(data_nf.x));
% 
% 
% surf(reshape(data_nf.x,numelX,numelY),...
%      reshape(data_nf.y,numelX,numelY),...
%      reshape(angle(data_nf.E(:,1))*180/pi,numelX,numelY))
% title(['Scan ' num2str(numelX) 'x' num2str(numelY) ' Phase X Component'])
% figure('name','Near-Field Phase','numbertitle','off',...
%         'units','normalized','outerposition',[0 0 1 1]);
% 
% numelX= sqrt(length(data_nf.x));
% numelY= sqrt(length(data_nf.x));
% 
% 
% surf(reshape(data_nf.x,numelX,numelY),...
%      reshape(data_nf.y,numelX,numelY),...
%      reshape(angle(data_nf.E(:,2))*180/pi,numelX,numelY))
% title(['Scan ' num2str(numelX) 'x' num2str(numelY) ' Phase Y Component'])
% figure('name','Near-Field Phase','numbertitle','off',...
%         'units','normalized','outerposition',[0 0 1 1]);
% 
% numelX= sqrt(length(data_nf.x));
% numelY= sqrt(length(data_nf.x));
% 
% 
% surf(reshape(data_nf.x,numelX,numelY),...
%      reshape(data_nf.y,numelX,numelY),...
%      reshape(angle(data_nf.E(:,3))*180/pi,numelX,numelY))
% title(['Scan ' num2str(numelX) 'x' num2str(numelY) ' Phase Z Component'])

end

