function [] = plotPlanarNFData(data_nf)

figure('name','Near-Field Magnitude','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1]);

numelX= sqrt(length(data_nf.x));
numelY= sqrt(length(data_nf.x));


surf(reshape(data_nf.x,numelX,numelY),...
     reshape(data_nf.y,numelX,numelY),...
     reshape(data_nf.Eabs,numelX,numelY))
title(['Scan ' num2str(numelX) 'x' num2str(numelY) 'Magnitude'])

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

