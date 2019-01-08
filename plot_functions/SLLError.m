function [sll_ff,sll_nf2ff,sll_err] = SLLError(data_ff,data_nf2ff,scanner_case)
%SLLERROR Side-Lobe-Level Error analysis


switch scanner_case
    case 'planar'
        phi_cut = pi/2;
        % Far-Field
        i = find(data_ff.phi==phi_cut);
        m = find(data_ff.phi==phi_cut+pi);
        ff_cut_angles = [-fliplr(data_ff.theta(m)') data_ff.theta(i)'];
        maxValue1 = (max(data_ff.Eabs(i)));
        maxValue2 = (max(data_ff.Eabs(m)));
        if maxValue1>maxValue2
            maxValue =maxValue1;
        else
            maxValue=maxValue2;
        end
        ff_cut = [fliplr(data_ff.Eabs(m)')'; data_ff.Eabs(i)]/maxValue;


        %  Near Field
        i = find(data_nf2ff.phi==phi_cut);
        nf2ff_cut_angles = data_nf2ff.theta(i)';
        maxValue = max(max(data_nf2ff.Eabs(i)));
        nf2ff_cut =  data_nf2ff.Eabs(i)/maxValue;
    
    case 'cylindrical'
        phi_cut=0;
        % Far-Field
        i = find(data_ff.phi==phi_cut);
        m = find(data_ff.phi==phi_cut+pi);
        ff_cut_angles = [-fliplr(data_ff.theta(m)') data_ff.theta(i)'];
        maxValue1 = (max(data_ff.Eabs(i)));
        maxValue2 = (max(data_ff.Eabs(m)));

        if maxValue1>maxValue2
            maxValue =maxValue1;
        else
            maxValue=maxValue2;
        end
        ff_cut = [fliplr(data_ff.Eabs(m)')'; data_ff.Eabs(i)]/maxValue;

        % Near-Field
        i = find(abs(data_nf2ff.phi-phi_cut)<0.000001);
        m = find(abs(data_nf2ff.phi-(phi_cut+pi))<0.000001);
        nf2ff_cut_angles = [-fliplr(data_nf2ff.theta(m)') data_nf2ff.theta(i)']; 
        maxValue1 = (max(data_nf2ff.Eabs(i)));
        maxValue2 = (max(data_nf2ff.Eabs(m)));

        if maxValue1>maxValue2
            maxValue =maxValue1;
        else
            maxValue=maxValue2;
        end
        nf2ff_cut = [fliplr(data_nf2ff.Eabs(m)')'; data_nf2ff.Eabs(i)]/maxValue;
    
    case 'spherical'
end


% Side-Lobe Level Far-Field
[pks_ff,idx_ff] = findpeaks(ff_cut);
[m,i] = max(pks_ff(pks_ff<0.5));
sll_ff = 20*log10(m);

if pks_ff(i) == 1
    i = i+1;
end
sll_nf2ff = 20*log10(nf2ff_cut(idx_ff(i)-90));

sll_err = sll_ff - sll_nf2ff;

end

