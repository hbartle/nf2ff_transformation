function [] = plotFFThetaCut(data_ff,theta,normalized,logarithmic)


i = find(data_ff.theta==theta);
ff_cut_angles = data_ff.phi(i); 
maxValue = max(data_ff.Eabs(i));
if normalized == true
    ff_cut = data_ff.Eabs(i)/maxValue;
else 
    ff_cut = data_ff.Eabs(i);
end

if logarithmic == true
    plot(ff_cut_angles*180/pi,20*log10(ff_cut))
elseif logarithmic == false
    plot(ff_cut_angles*180/pi,ff_cut)
end
hold on

end

