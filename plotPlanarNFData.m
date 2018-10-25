function [] = plotPlanarNFData(data_nf)

figure('name','Near-Field','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1]);

numelX= sqrt(length(data_nf.x));
numelY= sqrt(length(data_nf.x));


surf(reshape(data_nf.x,numelX,numelY),...
     reshape(data_nf.y,numelX,numelY),...
     reshape(data_nf.Eabs,numelX,numelY))
title(['Scan ' num2str(numelX) 'x' num2str(numelY)])


end

