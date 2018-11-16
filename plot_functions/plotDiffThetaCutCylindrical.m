function [] = plotDiffThetaCutCylindrical(data_nf2ff,data_ff,theta)
phi_range = unique(data_nf2ff.phi);


% Normalized Far-Field cut

i = find(data_ff.theta==theta & ismember(round(data_ff.phi*180/pi,0),round(phi_range*180/pi,0)));
ff_cut_angles = data_ff.phi(i); 
ff_cut =  data_ff.Eabs(i)/max(data_ff.Eabs);


% Normalized NF2FF cut
i = find(abs(data_nf2ff.theta-theta)<0.000001);

nf2ff_cut_angles = data_nf2ff.phi(i); 
maxValue = max(max(data_nf2ff.Eabs(i)));
nf2ff_cut = data_nf2ff.Eabs(i)/maxValue;


plot(nf2ff_cut_angles*180/pi,20*log10(abs(ff_cut - nf2ff_cut)))
hold on
xlim([0 360])
end

