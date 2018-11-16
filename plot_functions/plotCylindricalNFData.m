function [] = plotCylindricalNFData(data_nf)

figure('name','Near-Field Magnitude','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1]);

h = unique(data_nf.z);
phi = unique(round(atan2(data_nf.y,data_nf.x),4));
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

end

