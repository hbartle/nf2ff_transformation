function [] = plotFFPhiCut(data_ff,phi,normalized,logarithmic)


i = find(data_ff.phi==phi);
m = find(data_ff.phi==phi+pi);
ff_cut_angles = [-fliplr(data_ff.theta(m)') data_ff.theta(i)']; 
if normalized == true
    ff_cut = [fliplr(data_ff.Eabs(m)')'; data_ff.Eabs(i)]/max(data_ff.Eabs);
else 
    ff_cut = [fliplr(data_ff.Eabs(m)')'; data_ff.Eabs(i)];
end

if logarithmic == true
    plot(ff_cut_angles*180/pi,20*log10(ff_cut))
elseif logarithmic == false
    plot(ff_cut_angles*180/pi,ff_cut)
end
hold on

end

