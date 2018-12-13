function [] = plotDiffPhiCutCylindrical(data_nf2ff,data_ff,phi,theta_range)


% Normalized Far-Field cut
i = find(data_ff.phi==phi & ismember(data_ff.theta,theta_range));
m = find(data_ff.phi==phi+pi & ismember(data_ff.theta,theta_range));
ff_cut_angles = [-fliplr(data_ff.theta(m)') data_ff.theta(i)'];
maxValue1 = (max(data_ff.Eabs(i)));
maxValue2 = (max(data_ff.Eabs(m)));

if maxValue1>maxValue2
    maxValue =maxValue1;
else
    maxValue=maxValue2;
end

ff_cut = [fliplr(data_ff.Eabs(m)')'; data_ff.Eabs(i)]/maxValue;


% Normalized NF2FF cut
i = find(abs(data_nf2ff.phi-phi)<0.000001);
m = find(abs(data_nf2ff.phi-(phi+pi))<0.000001);
nf2ff_cut_angles = [-fliplr(data_nf2ff.theta(m)') data_nf2ff.theta(i)']; 
maxValue = max(max(data_nf2ff.Eabs(i)));
nf2ff_cut = [fliplr(data_nf2ff.Eabs(m)')'; data_nf2ff.Eabs(i)]/maxValue;


plot(nf2ff_cut_angles*180/pi,20*log10(abs(ff_cut - nf2ff_cut)))
hold on
end

