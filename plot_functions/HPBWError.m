function [hpbw_ff,hpbw_nf2ff,hpbw_err] = HPBWError(data_ff,data_nf2ff,scanner_case)
%HPBWERROR Calculate the error in the Half-Power Beamwidth of the antenna
%pattern.


%% Phi=0 cutting plane
phi_cut = 0;

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


% Near Field
switch scanner_case
    case 'planar'
        i = find(data_nf2ff.phi==phi_cut);
        nf2ff_cut_angles = data_nf2ff.theta(i)';
        maxValue = max(max(data_nf2ff.Eabs(i)));
        nf2ff_cut =  data_nf2ff.Eabs(i)/maxValue;

    case 'cylindrical'
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


% Find FF-HPBW 
i = find(ff_cut>0.5);
angle_res = mean(diff(ff_cut_angles));
hpbw_ff(1) = (ff_cut_angles(i(end)) - ff_cut_angles(i(1)))*180/pi;
hpbw_ff(1) = length(i)*angle_res*180/pi;

% Find NF-HPBW
i = find(nf2ff_cut>0.5);
angle_res = mean(diff(nf2ff_cut_angles));
hpbw_nf2ff(1) = (nf2ff_cut_angles(i(end)) - nf2ff_cut_angles(i(1)))*180/pi;
hpbw_nf2ff(1) = length(i)*angle_res*180/pi;


%% Phi=90 / Theta = 90

% Near Field
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
        
        % Near-Field
        i = find(data_nf2ff.phi==phi_cut);
        nf2ff_cut_angles = data_nf2ff.theta(i)';
        maxValue = max(max(data_nf2ff.Eabs(i)));
        nf2ff_cut =  data_nf2ff.Eabs(i)/maxValue;

        
        
    case 'cylindrical'
        % Cut in Theta=90 plane since antenna is directed in positive
        % x-axis
       
        % Far-Field
        theta_cut=pi/2;
        i = find(data_ff.theta==theta_cut);
        ff_cut_angles = data_ff.phi(i); 
        maxValue = max(data_ff.Eabs(i));
        ff_cut = data_ff.Eabs(i)/maxValue;
         

        % Near-field
        % Find values in cutting plane
        i = find(abs(data_nf2ff.theta-theta_cut)<0.000001);
        nf2ff_cut_angles = data_nf2ff.phi(i); 
        maxValue = max(max(data_nf2ff.Eabs(i)));
        nf2ff_cut = data_nf2ff.Eabs(i)/maxValue;
      
        
    case 'spherical'
end

% Find FF-HPBW 
i = find(ff_cut>0.5);
angle_res = mean(diff(ff_cut_angles));
% hpbw_ff(2) = (ff_cut_angles(i(end)) - ff_cut_angles(i(1)))*180/pi;
hpbw_ff(2) = length(i)*angle_res*180/pi;

% Find NF-HPBW
i = find(nf2ff_cut>0.5);
angle_res = mean(diff(nf2ff_cut_angles));
% hpbw_nf2ff(2) = (nf2ff_cut_angles(i(end)) - nf2ff_cut_angles(i(1)))*180/pi;
hpbw_nf2ff(2) = length(i)*angle_res*180/pi;


%% Find HPBW Error
hpbw_err = hpbw_ff - hpbw_nf2ff;

% Transpose results
hpbw_ff = hpbw_ff';
hpbw_nf2ff = hpbw_nf2ff';
hpbw_err = hpbw_err';
end

