function [] = plotNFPhiCutCylindrical(data_nf2ff,phi,normalized,logarithmic)

% Find values in cutting plane
i = find(abs(data_nf2ff.phi-phi)<0.000001);
m = find(abs(data_nf2ff.phi-(phi+pi))<0.000001);
nf2ff_cut_angles = [-fliplr(data_nf2ff.theta(m)') data_nf2ff.theta(i)']; 
maxValue1 = max(max(data_nf2ff.Eabs(i)));
maxValue2 = max(max(data_nf2ff.Eabs(m)));

if maxValue1>maxValue2
    maxValue =maxValue1;
else
    maxValue=maxValue2;
end

if normalized == true
    nf2ff_cut = [fliplr(data_nf2ff.Eabs(m)')'; data_nf2ff.Eabs(i)]/maxValue;
else 
    nf2ff_cut = [fliplr(data_nf2ff.Eabs(m)')'; data_nf2ff.Eabs(i)];
end

if logarithmic == true
    plot(nf2ff_cut_angles*180/pi,20*log10(nf2ff_cut))
elseif logarithmic == false
    plot(nf2ff_cut_angles*180/pi,nf2ff_cut)
end
hold on
end

