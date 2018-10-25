function [] = plotDiffPhiCut(data_nf2ff,data_ff,phi,theta_range)


% Normalized Far-Field cut
i = find(data_ff.phi==phi & ismember(data_ff.theta,theta_range));
m = find(data_ff.phi==phi+180 & ismember(data_ff.theta,theta_range));
ff_cut_angles = [-fliplr(data_ff.theta(m)') data_ff.theta(i)']; 
ff_cut = [fliplr(data_ff.Eabs(m)')'; data_ff.Eabs(i)]/max(data_ff.Eabs);


% Normalized NF2FF cut
i = find(data_nf2ff.phi==phi);
m = find(data_nf2ff.phi==phi+180);

nf2ff_cut_angles = [-fliplr(data_nf2ff.theta(m)') data_nf2ff.theta(i)'];
maxValue = max(max(data_nf2ff.Eabs(m)),max(data_nf2ff.Eabs(i)));
nf2ff_cut = [fliplr(data_nf2ff.Eabs(m)')'; data_nf2ff.Eabs(i)]/maxValue;

plot(nf2ff_cut_angles,20*log10(abs(ff_cut - nf2ff_cut)))
hold on
end

