function [] = plotNFPhiCut(data_nf2ff,phi,normalized)


% Find values in cutting plane
i = find(data_nf2ff.phi==phi);
m = find(data_nf2ff.phi==phi+180);

% Flip one have of the values to make continous plot
nf2ff_cut_angles = [-fliplr(data_nf2ff.theta(m)') data_nf2ff.theta(i)'];
maxValue = max(max(data_nf2ff.Eabs(m)),max(data_nf2ff.Eabs(i)));
if normalized == true
    nf2ff_cut = [fliplr(data_nf2ff.Eabs(m)')'; data_nf2ff.Eabs(i)]/maxValue;
else
    nf2ff_cut = [fliplr(data_nf2ff.Eabs(m)')'; data_nf2ff.Eabs(i)];
end
plot(nf2ff_cut_angles,nf2ff_cut)

end

